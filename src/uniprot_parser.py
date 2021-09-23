import requests
from tqdm import tqdm
from pathlib import Path

# Cross References
DB_ABBREVS = ["database(EMBL)"] + ["database(" + line.strip().split(": ")[1] + ")" for line in
                                   requests.get("https://www.uniprot.org/docs/dbxref.txt").text.split("\n") if
                                   "Abbrev:" in line]

COLUMN_NAMES = [
    # Names & Taxonomy
    'id', 'entry name', 'genes', 'genes(PREFERRED)', 'genes(ALTERNATIVE)',
    'genes(OLN)', 'genes(ORF)', 'organism', 'organism-id', 'protein names',
    'proteome', 'lineage(ALL)', 'lineage-id', 'virus hosts',
    # Sequences
    'fragement', 'sequence', 'length', 'mass', 'encodedon',
    'comment(ALTERNATIVE PRODUCTS)', 'comment(ERRONEOUS GENE MODEL PREDICTION)',
    'comment(ERRONEOUS INITIATION)', 'comment(ERRONEOUS TERMINATION)',
    'comment(ERRONEOUS TRANSLATION)', 'comment(FRAMESHIFT)',
    'comment(MASS SPECTROMETRY)', 'comment(POLYMORPHISM)',
    'comment(RNA EDITING)', 'comment(SEQUENCE CAUTION)',
    'feature(ALTERNATIVE SEQUENCE)', 'feature(NATURAL VARIANT)',
    'feature(NON ADJACENT RESIDUES)',
    'feature(NON STANDARD RESIDUE)', 'feature(NON TERMINAL RESIDUE)',
    'feature(SEQUENCE CONFLICT)', 'feature(SEQUENCE UNCERTAINTY)',
    'version(sequence)',
    # Family and Domains
    'domains', 'domain', 'comment(DOMAIN)', 'comment(SIMILARITY)',
    'feature(COILED COIL)', 'feature(COMPOSITIONAL BIAS)',
    'feature(DOMAIN EXTENT)', 'feature(MOTIF)', 'feature(REGION)',
    'feature(REPEAT)', 'feature(ZINC FINGER)',
    # Function
    'ec', 'comment(ABSORPTION)', 'comment(CATALYTIC ACTIVITY)',
    'comment(COFACTOR)', 'comment(ENZYME REGULATION)', 'comment(FUNCTION)',
    'comment(KINETICS)', 'comment(PATHWAY)', 'comment(REDOX POTENTIAL)',
    'comment(TEMPERATURE DEPENDENCE)', 'comment(PH DEPENDENCE)',
    'feature(ACTIVE SITE)', 'feature(BINDING SITE)', 'feature(DNA BINDING)',
    'feature(METAL BINDING)', 'feature(NP BIND)', 'feature(SITE)',
    # Gene Ontologys
    'go', 'go(biological process)', 'go(molecular function)',
    'go(cellular component)', 'go-id',
    # InterPro
    'interpro',
    # Interaction
    'interactor', 'comment(SUBUNIT)',
    # Publications
    'citation', 'citationmapping',
    # Date of
    'created', 'last-modified', 'sequence-modified', 'version(entry)',
    # Structure
    '3d', 'feature(BETA STRAND)', 'feature(HELIX)', 'feature(TURN)',
    # Subcellular location
    'comment(SUBCELLULAR LOCATION)', 'feature(INTRAMEMBRANE)',
    'feature(TOPOLOGICAL DOMAIN)',
    'feature(TRANSMEMBRANE)',
    # Miscellaneous
    'annotation score', 'score', 'features', 'comment(CAUTION)',
    'comment(TISSUE SPECIFICITY)',
    'comment(GENERAL)', 'keywords', 'context', 'existence', 'tools',
    'reviewed', 'feature', 'families', 'subcellular locations', 'taxonomy',
    'version', 'clusters', 'comments', 'database', 'keyword-id', 'pathway',
    'score',
    # Pathology & Biotech
    'comment(ALLERGEN)', 'comment(BIOTECHNOLOGY)', 'comment(DISRUPTION PHENOTYPE)',
    'comment(DISEASE)', 'comment(PHARMACEUTICAL)', 'comment(TOXIC DOSE)',
    # PTM / Processsing
    'comment(PTM)', 'feature(CHAIN)', 'feature(CROSS LINK)', 'feature(DISULFIDE BOND)',
    'feature(GLYCOSYLATION)', 'feature(INITIATOR METHIONINE)', 'feature(LIPIDATION)',
    'feature(MODIFIED RESIDUE)', 'feature(PEPTIDE)', 'feature(PROPEPTIDE)',
    'feature(SIGNAL)', 'feature(TRANSIT)',
    # Taxonomic lineage
    'lineage(all)', 'lineage(SUPERKINGDOM)', 'lineage(KINGDOM)', 'lineage(SUBKINGDOM)',
    'lineage(SUPERPHYLUM)', 'lineage(PHYLUM)', 'lineage(SUBPHYLUM)', 'lineage(SUPERCLASS)',
    'lineage(CLASS)', 'lineage(SUBCLASS)', 'lineage(INFRACLASS)', 'lineage(SUPERORDER)',
    'lineage(ORDER)', 'lineage(SUBORDER)', 'lineage(INFRAORDER)', 'lineage(PARVORDER)',
    'lineage(SUPERFAMILY)', 'lineage(FAMILY)', 'lineage(SUBFAMILY)', 'lineage(TRIBE)',
    'lineage(SUBTRIBE)', 'lineage(GENUS)', 'lineage(SUBGENUS)', 'lineage(SPECIES GROUP)',
    'lineage(SPECIES SUBGROUP)', 'lineage(SPECIES)', 'lineage(SUBSPECIES)', 'lineage(VARIETAS)',
    'lineage(FORMA)',
    # Taxonomic identifier
    'lineage-id(all)', 'lineage-id(SUPERKINGDOM)', 'lineage-id(KINGDOM)', 'lineage-id(SUBKINGDOM)',
    'lineage-id(SUPERPHYLUM)', 'lineage-id(PHYLUM)', 'lineage-id(SUBPHYLUM)', 'lineage-id(SUPERCLASS)',
    'lineage-id(CLASS)', 'lineage-id(SUBCLASS)', 'lineage-id(INFRACLASS)', 'lineage-id(SUPERORDER)',
    'lineage-id(ORDER)', 'lineage-id(SUBORDER)', 'lineage-id(INFRAORDER)', 'lineage-id(PARVORDER)',
    'lineage-id(SUPERFAMILY)', 'lineage-id(FAMILY)', 'lineage-id(SUBFAMILY)', 'lineage-id(TRIBE)',
    'lineage-id(SUBTRIBE)', 'lineage-id(GENUS)', 'lineage-id(SUBGENUS)', 'lineage-id(SPECIES GROUP)',
    'lineage-id(SPECIES SUBGROUP)', 'lineage-id(SPECIES)', 'lineage-id(SUBSPECIES)', 'lineage-id(VARIETAS)',
    'lineage-id(FORMA)']


def get_uniprot_info_from_ids(ids: list, filename: Path, chunk=False, identifier: str = "ACC+ID", to: str = "ACC",
                              columns: str = ",".join(COLUMN_NAMES)):
    """
    Batch retrieval of IDs and information from UniProt.

    Parameters
    ----------
    ids
        input IDs
    filename
        write to this file
    chunk
        split into multiple queries of size 100 and join results
    identifier
        type of input IDs
    to
        output ID format - ACC returns all column information
    columns
        column names to return, preformatted sting (",".join(column_names))


    Returns
    -------
    All information written as newline separated, tab delimited text.
    """
    mapping_url = 'http://www.uniprot.org/uploadlists/'
    mapping_params = {
        'from': identifier,
        'to': to,
        'format': 'tab',
        'columns': columns
    }
    with open(filename, "w") as f:
        if chunk:
            num_tries = 5
            for i, id_i in tqdm(enumerate(range(0, len(ids), 100))):
                id_chunk = ids[id_i: id_i + 100]
                good_text = False
                text = ""
                try_number = 0
                while not good_text and try_number < num_tries:
                    mapping_params['query'] = ' '.join(id_chunk)
                    response = requests.post(mapping_url, params=mapping_params)
                    text = response.text
                    if "<html><head>" in text:
                        good_text = False
                    else:
                        good_text = True
                    try_number += 1
                if i == 0:
                    f.write(text)
                else:
                    f.write("\n".join(text.split("\n")[1:]))
        else:
            mapping_params['query'] = ' '.join(ids)
            response = requests.post(mapping_url, params=mapping_params)
            f.write(response.text)
