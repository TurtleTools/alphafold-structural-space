from pathlib import Path
from ftplib import FTP
import prody as pd
from dataclasses import dataclass
import numpy as np
import typing as ty
from geometricus import MomentInvariants, SplitType
import tarfile
from time import time
from tqdm import tqdm
from scipy import ndimage

from src import uniprot_parser, proteinnet_parser

UNIPROT_COLUMNS = ",".join(("id", "entry name", 'genes', 'genes(PREFERRED)', 'genes(ALTERNATIVE)',
                        'genes(OLN)', 'genes(ORF)', "organism", "protein names", "families",
                        'go', 'go(biological process)', 'go(molecular function)',
                        'go(cellular component)', 'database(PDB)', 'database(Pfam)'))


@dataclass
class MomentInvariantsSavable:
    name: str
    moments: ty.Union[np.ndarray, None]
    coordinates: ty.Union[np.ndarray, None]

    @classmethod
    def from_invariant(cls, invariant: MomentInvariants):
        return cls(invariant.name, invariant.moments, invariant.coordinates)


def download_data(output_folder: Path):
    if not output_folder.exists():
        output_folder.mkdir()
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    ftp.cwd("/pub/databases/alphafold")
    for filename in ftp.nlst():
        print(f"Retrieving {filename}")
        with open(output_folder / filename, 'wb') as f:
            ftp.retrbinary('RETR ' + filename, f.write)


def extract_data(input_folder, output_folder):
    for filename in input_folder.glob("*.tar"):
        tar = tarfile.open(str(filename))
        folder = output_folder / str(filename.stem)
        if not folder.exists():
            folder.mkdir()
        tar.extractall(str(folder))  # specify which folder to extract to
        tar.close()


def get_uniprot_info(data_folder, aux_folder, extension="pdb.gz"):
    start_time = time()
    for folder in data_folder.iterdir():
        if folder.is_dir() and folder.stem.startswith("UP0"):
            uniprot_file = aux_folder / f"{folder.stem}_uniprot.txt"
            uniprot_ids = [filename.stem.split("-")[1] for filename in folder.glob(f"*{extension}")]
            if not uniprot_file.exists():
                uniprot_parser.get_uniprot_info_from_ids(uniprot_ids, uniprot_file, chunk=True,
                                                         columns=UNIPROT_COLUMNS)
            print(f"{folder.stem}: Time elapsed: {time() - start_time}s")


def get_AF_shapemers(root_folder,
                     resolution_kmer=4,
                     resolution_radius=6,
                     length_threshold=50):
    root_folder = Path(root_folder)
    with open(
            root_folder /
            f"AF_ids_corpus_resolution_{resolution_kmer}_{resolution_radius}_threshold_{length_threshold}.txt",
            "w") as corpus_file:
        for folder in root_folder.iterdir():
            if folder.is_dir() and folder.stem.startswith("UP0"):
                print(folder.stem)
                for i, pdb_file in tqdm(enumerate(folder.glob("*.pdb.gz"))):
                    key = pdb_file.stem
                    pdb = pd.parsePDB(str(pdb_file)).select("protein and calpha")
                    betas = ndimage.gaussian_filter1d(pdb.getBetas(), sigma=5)
                    coords = pdb.getCoords()
                    sequence = pdb.getSequence()

                    indices = np.ones(betas.shape[0], dtype=int)
                    indices[np.where(betas < 70)] = 0

                    slices = ndimage.find_objects(ndimage.label(indices)[0])
                    index = 0
                    shapemers = []
                    for s in slices:
                        s = s[0]
                        if s.stop - s.start > length_threshold:
                            index += 1
                            invariants = MomentInvariants.from_coordinates(
                                key,
                                coords[s.start: s.stop],
                                sequence[s.start: s.stop],
                                split_type=SplitType.KMER_CUT,
                                split_size=16
                            )
                            shapemers += [f"k{x[0]}i{x[1]}i{x[2]}i{x[3]}" for x in
                                          (np.log1p(invariants.moments) * resolution_kmer).astype(int)]
                            invariants = MomentInvariants.from_coordinates(
                                key,
                                coords[s.start: s.stop],
                                sequence[s.start: s.stop],
                                split_type=SplitType.RADIUS,
                                split_size=10
                            )
                            shapemers += [f"r{x[0]}i{x[1]}i{x[2]}i{x[3]}" for x in
                                          (np.log1p(invariants.moments) * resolution_radius).astype(int)]
                    if index > 0:
                        corpus_file.write(key + "\t" + " ".join(shapemers) + "\n")


def get_PDB_shapemers(casp_file, root_folder, resolution_kmer=4, resolution_radius=6):
    with open(Path(root_folder) / f"PDB_{casp_file.stem}_ids_corpus_resolution_{resolution_kmer}_{resolution_radius}.txt",
              "w") as corpus_file:
        for entry in tqdm(proteinnet_parser.yield_records_from_file(casp_file, 20)):
            entry = proteinnet_parser.clean_entry(entry, 'ca')
            invariants = MomentInvariants.from_coordinates(
                entry["ID"],
                entry["tertiary"],
                entry["primary"],
                split_type=SplitType.KMER_CUT,
                split_size=16
            )
            shapemers = [f"k{x[0]}i{x[1]}i{x[2]}i{x[3]}" for x in
                         (np.log1p(invariants.moments) * resolution_kmer).astype(int)]
            invariants = MomentInvariants.from_coordinates(
                entry["ID"],
                entry["tertiary"],
                entry["primary"],
                split_type=SplitType.RADIUS,
                split_size=10
            )
            shapemers += [f"r{x[0]}i{x[1]}i{x[2]}i{x[3]}" for x in
                          (np.log1p(invariants.moments) * resolution_radius).astype(int)]
            if len(shapemers):
                corpus_file.write(entry["ID"] + "\t" + " ".join(shapemers) + "\n")


def get_AF_protein_information(data_folder):
    data_folder = Path(data_folder)
    avg_scores = {}
    lengths_high_confidence = {}
    lengths_full = {}
    for folder in data_folder.iterdir():
        if folder.is_dir() and folder.stem.startswith("UP0"):
            for filename in tqdm(folder.glob("*.pdb.gz")):
                pdb = pd.parsePDB(str(filename))
                if pdb is None:
                    continue
                key = filename.stem
                avg_scores[key] = np.median(pdb.getBetas())
                pdb = pdb.select("protein and calpha")
                if pdb is None:
                    continue
                lengths_full[key] = len(pdb)
                pdb = pdb.select("beta > 70")
                if pdb is None:
                    continue
                lengths_high_confidence[key] = len(pdb)
    return avg_scores, lengths_high_confidence, lengths_full


def get_PDB_protein_information(casp_files):
    coords = {}
    for casp_file in casp_files:
        for entry in tqdm(proteinnet_parser.yield_records_from_file(casp_file, 20)):
            entry = proteinnet_parser.clean_entry(entry, 'ca')
            coords[entry["ID"]] = entry["tertiary"]
    return coords
