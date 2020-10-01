import os
import py2neo
import logging
import json

logging.basicConfig(level=logging.INFO)
logging.getLogger('biomedgraph').setLevel(level=logging.DEBUG)

from biomedgraph.datasources import Gtex, GeneOntology, Reactome, NcbiTaxonomy, NcbiGene, BigWordList, Ensembl, Refseq, \
    Uniprot
from biomedgraph.parser import *

log = logging.getLogger(__name__)

ROOT_DIR = os.getenv('ROOT_DIR', '/download')
RUN_MODE = os.getenv('RUN_MODE', 'prod')

NEO4J_CONFIG_STRING = os.getenv("NEO4J")
log.info(NEO4J_CONFIG_STRING)

try:
    NEO4J_CONFIG_DICT = json.loads(NEO4J_CONFIG_STRING)
except json.decoder.JSONDecodeError:

    NEO4J_CONFIG_STRING = NEO4J_CONFIG_STRING.replace("'", '"')
    log.info(NEO4J_CONFIG_STRING)
    NEO4J_CONFIG_DICT = json.loads(NEO4J_CONFIG_STRING)

log.info(NEO4J_CONFIG_DICT)

def run_parser(parser):
    """
    Run a parser and log.

    :param parser: The parser
    :return: The parser after running
    """
    log.info("Run parser {}".format(parser.__class__.__name__))
    parser.run_with_mounted_arguments()
    log.info(parser.container.nodesets)
    log.info(parser.container.relationshipsets)
    return parser


def create_index(graph, parser_list):
    """
    Create all indexes for the RelationshipSets in a list of Parsers.

    :param graph: Py2neo graph instance
    :param parser_list: List of parsers
    """
    for parser in parser_list:
        for relationshipset in parser.container.relationshipsets:
            relationshipset.create_index(graph)
        for nodeset in parser.container.nodesets:
            nodeset.create_index(graph)


def create_nodesets(graph, parser_list):
    """
    Create the NodeSets for a list of parsers

    :graph: Py2neo graph instance
    :param parser_list: List of Parsers
    """
    for parser in parser_list:
        log.info("Create nodes for parser {}".format(parser.__class__.__name__))
        for nodeset in parser.container.nodesets:
            nodeset.merge(graph)


def create_relationshipsets(graph, parser_list):
    """
    Create the RelationshipSets for a list of parsers

    :graph: Py2neo graph instance
    :param parser_list: List of Parsers
    """
    for parser in parser_list:
        log.info("Create relationships for parser {}".format(parser.__class__.__name__))
        for relset in parser.container.relationshipsets:
            relset.merge(graph)


if __name__ == '__main__':

    if RUN_MODE.lower() == 'test':
        log.info("Run tests")

    else:
        graph = py2neo.Graph(host=NEO4J_CONFIG_DICT['host'], user=NEO4J_CONFIG_DICT['user'], password=NEO4J_CONFIG_DICT['password'], secure=NEO4J_CONFIG_DICT['secure'], verify=False)

        # Download Datasources
        # ====================
        log.info('Download NCBI Gene')
        ncbigene = NcbiGene(ROOT_DIR)
        if not ncbigene.latest_local_instance():
            ncbigene.download()

        log.info('Download Big Word List')
        bigwordlist = BigWordList(ROOT_DIR)
        if not bigwordlist.latest_local_instance():
            bigwordlist.download()

        log.info('Download ENSEMBL')
        ensembl = Ensembl(ROOT_DIR)
        if not ensembl.latest_local_instance():
            ensembl.download(ensembl.latest_remote_version(), taxids=['9606'])

        log.info('Download RefSeq')
        refseq = Refseq(ROOT_DIR)
        if not refseq.latest_local_instance():
            refseq.download(refseq.latest_remote_version())

        log.info('Download Uniprot')
        uniprot = Uniprot(ROOT_DIR)
        if not uniprot.latest_local_instance():
            uniprot.download()

        gtex = Gtex(ROOT_DIR)
        if not gtex.latest_local_instance():
            gtex.download()

        reactome = Reactome(ROOT_DIR)
        if not reactome.latest_local_instance():
            reactome.download()

        ncbi_taxonomy = NcbiTaxonomy(ROOT_DIR)
        if not ncbi_taxonomy.latest_local_instance():
            ncbi_taxonomy.download()

        geneontology = GeneOntology(ROOT_DIR)
        if not geneontology.latest_local_instance():
            geneontology.download(taxids=['9606'])

        # run Parsers
        # ================
        parsers_done = []

        ncbigene_parser = NcbiGeneParser(ROOT_DIR)
        ncbigene_parser.datasource_instances.append(ncbigene.latest_local_instance())
        ncbigene_parser.taxid = '9606'
        parsers_done.append(run_parser(ncbigene_parser))

        bigwordlist_parser = BigWordListParser(ROOT_DIR)
        bigwordlist_parser.datasource_instances.append(bigwordlist.latest_local_instance())
        parsers_done.append(run_parser(bigwordlist_parser))

        ensembl_entity_parser = EnsemblEntityParser(ROOT_DIR)
        ensembl_entity_parser.datasource_instances.append(ensembl.latest_local_instance())
        ensembl_entity_parser.taxid = '9606'
        parsers_done.append(run_parser(ensembl_entity_parser))

        ensembl_mapping_parser = EnsemblMappingParser(ROOT_DIR)
        ensembl_mapping_parser.datasource_instances.append(ensembl.latest_local_instance())
        ensembl_mapping_parser.taxid = '9606'
        parsers_done.append(run_parser(ensembl_mapping_parser))

        refseq_entity_parser = RefseqEntityParser(ROOT_DIR)
        refseq_entity_parser.datasource_instances.append(refseq.latest_local_instance())
        refseq_entity_parser.taxid = '9606'
        parsers_done.append(run_parser(refseq_entity_parser))

        refseq_codes_parser = RefseqCodesParser(ROOT_DIR)
        refseq_codes_parser.datasource_instances.append(refseq.latest_local_instance())
        refseq_codes_parser.taxid = '9606'
        parsers_done.append(run_parser(refseq_codes_parser))

        uniprot_knowledgebase_parser = UniprotKnowledgebaseParser(ROOT_DIR)
        uniprot_knowledgebase_parser.datasource_instances.append(uniprot.latest_local_instance())
        uniprot_knowledgebase_parser.taxid = '9606'
        parsers_done.append(run_parser(uniprot_knowledgebase_parser))

        go_annotation_parser = GeneOntologyAssociationParser(ROOT_DIR)
        go_annotation_parser.datasource_instances.append(geneontology.latest_local_instance())
        go_annotation_parser.taxid = '9606'
        parsers_done.append(run_parser(go_annotation_parser))

        gtex_parser = GtexMetadataParser(ROOT_DIR)
        gtex_parser.datasource_instances.append(gtex.latest_local_instance())
        parsers_done.append(run_parser(gtex_parser))

        gtex_data_parser = GtexDataParser(ROOT_DIR)
        gtex_data_parser.datasource_instances.append(gtex.latest_local_instance())
        parsers_done.append(run_parser(gtex_data_parser))

        reactome_parser = ReactomePathwayParser(ROOT_DIR)
        reactome_parser.datasource_instances.append(reactome.latest_local_instance())
        reactome_parser.datasource_instances.append(ncbi_taxonomy.latest_local_instance())
        parsers_done.append(run_parser(reactome_parser))

        reactome_mapping_parser = ReactomeMappingParser(ROOT_DIR)
        reactome_mapping_parser.datasource_instances.append(reactome.latest_local_instance())
        reactome_mapping_parser.datasource_instances.append(ncbi_taxonomy.latest_local_instance())
        reactome_mapping_parser.taxid = '9606'
        parsers_done.append(run_parser(reactome_mapping_parser))

        # Load data
        # ================
        create_index(graph, parsers_done)
        create_nodesets(graph, parsers_done)
        create_relationshipsets(graph, parsers_done)
