#!/usr/bin/env python3
import random
import click
#import utils as CMD
import not_alike3.utils as CMD
#import biostr as BS
import not_alike3.biostr as BS

@click.group()
def main():
    """
        Not-Alike: Command pipeline that identifies dissimilar regions of a target genome by comparing it to a genomes database.
    """
    pass

######################################
######          SEARCH          ######
######################################

@main.command()
@click.option('-g', '--genome', \
                help = 'Query genome FASTA file name', \
                required = True, \
                type = str)
@click.option('-db', '--database-file', \
                help = 'A name of the text file which contains the name of BLAST (*.db) files.', \
                required = True, \
                type = str)
@click.option('-c', '--comment', \
                help = 'Leave a comment enclosed by single quotes', \
                required = False, \
                default = 'Empty', \
                type = str)
@click.option('--config-file', \
                help = 'TOML Configuration file.', \
                required = True, \
                type = str)
def search(genome, database_file, comment, config_file):
    """
        Searches for not alike fragments in query genome
    """
    CMD.check_path_exists('split_out')
    CMD.check_path_exists('blast_db')
    CMD.check_path_exists('blast_out')
    CMD.check_path_exists('ht2_idx')
    CMD.check_path_exists('mapping')
    CMD.check_path_exists('gtfs')
    CMD.check_path_exists('log')

    # PARSING TOML
    args = CMD.parsing_toml(config_file)


    window_size = args['split'][0]['window_size']
    step_size = args['split'][0]['step_size']
    num_cores = args['global']['num_of_cores']
#    quite_opposite = args['global']['quite_opposite']

    PID = random.randrange(1, 9999999999)
    out_split = '.'.join(genome.split('.')[:-1]) + '_split_' + str(window_size) + '_' + str(step_size) + '.fasta'
    out_split = ''.join(out_split.split('/')[-1])
    out_split = 'split_out/' + out_split

    input_split = 'input_split.' + str(PID) + '.fasta'
    input_split = 'split_out/' + input_split

    print(out_split + ' was loaded.')
    print(input_split + ' was loaded.')

    if not CMD.os.path.exists(out_split):
        '''
            loading and splitting genome
        '''
        seqs = BS.loadDNASeqs(genome)
        sptSeqs = BS.splitBioString(seqs, window_size, step_size)
        seqs.freeDNA()
        sptSeqs.writeNoHideToFile(out_split)
    else:
        print('A split-genome file was found!!!')
        #B80. It didn't load DNA Seqs when it founds a split sequences file.
        # Bug solved.
        sptSeqs = BS.loadDNASeqs(out_split)

    CMD.copy_file(out_split, input_split)

    database_file_path = '/'.join(database_file.split('/')[:-1])

    db_files = CMD.load_lines(database_file)

    for search in args['search']:
        identity = search['identity']
        qcov = search['qcov']
        evalue = search['evalue']
        task = search['task']
        name = search['name']
        quite_opposite = search['quite_opposite']
    
        CMD.loggin_data(PID, genome, \
                        database_file, window_size, \
                        step_size, config_file, \
                        identity, qcov, \
                        evalue, task, \
                        name, comment, \
                        num_cores, quite_opposite)

        for f in db_files:
            dbf = '.'.join(f.split('.')[:-1]) + '.db'
            print('Blasting ' + dbf + ' ...')
            
            CMD.do_blast(input_split, \
                            database_file_path + '/' + dbf, \
                            'blast_out/out.blast', \
                            evalue, identity, \
                            qcov, task, num_cores)

            print('Updating input_split')
#            print('Loading blast output qseqids')
            lkdl_headers = BS.loadLkdList('blast_out/out.blast')
#            print('Filtering Bioseqs')
            sptSeqs.filterBioseq(lkdl_headers)
#            print('Writting non-hidden sequences')
            sptSeqs.writeNoHideToFile(input_split)

        print(f'{name} search strategy finished')
        print('BLASTn searching done!')
        print('Mapping on process.')
    
        CMD.mapping(genome, input_split)

        print('Assembly on process.')
        CMD.assembly(PID)

        print('Extracting sequences.')
        CMD.extseq(genome, PID)

        print('Doing assembly stats')
        CMD.do_assembly_stats('gtfs/nal_frags.' + str(PID) + '.fasta', PID)

        
        PID = random.randrange(1, 9999999999)
        CMD.copy_file(input_split, 'split_out/tmp.split')

        input_split = 'input_split.' + str(PID) + '.fasta'
        input_split = 'split_out/' + input_split

        CMD.copy_file('split_out/tmp.split', input_split)

    sptSeqs.freeDNA()
    lkdl_headers.freeLkdList()

    print('not-alike has finished.')

##############################################
######          DB MAKEBLAST            ######
##############################################

@main.command()
@click.option('-db', '--db-path', \
                help = 'Path to FASTA files database', \
                required = True, \
                type = str)
def db_makeblast(db_path):
    """
        Builds a BLAST_DB (version 5) database files
    """

    CMD.make_db(db_path)
    
##########################################
######          DB MAKEFILE         ######
##########################################

@main.command()
@click.option('-qg', '--query-genome', \
                help = 'Path to query genome (FASTA)', \
                required = True, \
                type = str)
@click.option('-db', '--db-path', \
                help = 'Path to FASTA files database', \
                required = True, \
                type = str)
@click.option('-e', '--exclude', \
                help = 'A list of accession numbers from the organisms you want to exclude from database text file.', \
                required = False, \
                default = None, \
                type = str)
@click.option('-i', '--include', \
                help = 'A list of accession numbers from the organism you want to include in database text file.', \
                required = False, \
                default = None, \
                type = str)
@click.option('-o', '--out', \
                help = 'Output file name', \
                required = True, \
                type = str)
def db_makefile(query_genome, db_path, exclude, include, out):
    """
        Creates the database text file which contains the BLAST_DB files paths.
    """
    if exclude == None and include == None:
        print('ERROR, Exclude and include options are empty')
        return 1
    if exclude != None and include != None:
        print('ERROR, Only one option between exclude or include is allowed')
        return 1
    if exclude != None:
        exclude = exclude.split(',')
    if include != None:
        include = include.split(',')
    CMD.make_txtfiledb(query_genome, db_path, exclude, include, out)

######################################
######          SHOW DB         ######
######################################

@main.command()
@click.option('-p', '--db-path', \
                help = 'Genomes database path', \
                required = True, \
                type = str)
def show_db(db_path):
    """
        Shows metadata of database [accession number, organism name and organism taxon id]
    """

    assembly_report = db_path + '/assembly_data_report.jsonl'
    assembly_report_tsv = db_path + '/assembly_data_report.tsv'
    if not CMD.file_exists(assembly_report_tsv):
        CMD.dataformat_tsv(assembly_report, assembly_report_tsv)

    print(assembly_report_tsv)

    CMD.print_table(assembly_report_tsv)

##########################################
######          SHOW EXP            ######
##########################################

@main.command()
@click.option('--sort-by', \
                help = 'Criteria to sort values.', \
                required = False, \
                default = 'DATE', \
                type = str)
def show_exp(sort_by):
    """
        Shows information about epxeriments stored in the current working directory
    """
    CMD.show_exp_info(sort_by)

##########################################
######          ASSM STATS          ######
##########################################

@main.command()
@click.option('-f', '--file-name', \
                help = 'FASTA file name', \
                required = True, \
                type = str)
@click.option('-pid', \
                help = 'Process ID', \
                required = False, \
                default = 00, \
                type = int)
def assm_stats(file_name, pid):
    """
        Calculates assembly statistics such as: Mean, Median, N50 and L50.
    """
    CMD.do_assembly_stats(file_name, pid)

##############################################
######          PRIMER SELECT           ######
##############################################

@main.command()
@click.option('--input-file', \
                help = 'Input fasta file', \
                required = True, \
                type = str)
@click.option('--opt-size', \
                help = 'Optimum primer size (nt)', \
                required = True, \
                type = int)
@click.option('--opt-gc', \
                help = 'Optimum GC percentage (%)', \
                required = True, \
                type = float)
@click.option('--opt-tm', \
                help = 'Optimum melting point (°C)', \
                required = True, \
                type = float)
@click.option('--product-size', \
                help = 'Expected product zise (bp) [i.e. 75-100]', \
                required = True, \
                type = str)
@click.option('--template-size-range', \
                help = 'Template sequence size range (bp) [i.e. 750-1000]', \
                required = True, \
                type = str)
def search_primers(input_file, opt_size, opt_gc, opt_tm, product_size, template_size_range):
    '''
        Selects the best fitted primer sequences based on input arguments.
    '''
    CMD.find_primers(input_file, opt_size, opt_gc, opt_tm, product_size, template_size_range)

if __name__ == '__main__':
    main()
