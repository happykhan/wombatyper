"""
Script to fetch short reads (with known mlst) and check they give the same ST call
"""
import csv 
import re
import logging 
import os
import subprocess
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s')

def main(ref_file='test_bed/SarAB.csv'): 
    # Open data table, 
    records = {} 
    with open(ref_file) as f:
        for row in csv.DictReader(f, dialect=csv.excel_tab):
            row_name = re.split(r'\(|;|\=|\/', row['Name'])[0].strip()
            accession = row['Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Bases;Average Length;Status)'].split(';')[0]
            if not accession.startswith('traces'):
                records[row_name] = { 
                    "accession" : accession,
                    "ST": row['ST']
                }
                for col in ["aroC", "dnaN", "hemD", "hisD", "purE","sucA", "thrA"  ]:
                    records[row_name][col] = row[col]

        logging.info('Done reading file ')

        # For each row, Fetch reads 
        if not os.path.exists('reads'):
            os.mkdir('reads')
        count = 0 
        for name, record in records.items():
            consensus = f"{name}.fasta"
            acc = record['accession']
            output = f'reads/{acc}_1.fastq'
            r2_output = f'reads/{acc}_2.fastq' 
            bam_file = f"{name}.bam"                       
            if not os.path.exists(consensus):
                if not os.path.exists(output):
                    cmd = f'./sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump {acc} -o reads/{acc}'
                    logging.info('Running ' + cmd )
                    subprocess.call(cmd, shell=True) 
                # Map consensus to reference 
                ref = "test_bed/AE006468.2.fna"
                if not os.path.exists(bam_file):
                    cmd = f"minimap2 -ax sr {ref}   {output}  {r2_output} | samtools sort > {bam_file}"            
                    logging.info('Running ' + cmd )
                    subprocess.call(cmd, shell=True) 
                if not os.path.exists(consensus):
                    cmd = f"/usr/users/QIB_fr005/alikhan/miniconda3/envs/samtools/bin/samtools consensus {bam_file} > {consensus}" 
                    logging.info('Running ' + cmd )
                    subprocess.call(cmd, shell=True)    
            if os.path.exists(output):
                os.remove(output)
            if os.path.exists(r2_output):
                os.remove(r2_output)    
            if os.path.exists(bam_file):
                os.remove(bam_file)                                    
            mlst_file = f"{name}.mlst"            
            if not os.path.exists(mlst_file):
                cmd = f"mlst {name}.fasta > {mlst_file}" 
                logging.info('Running ' + cmd )
                subprocess.call(cmd, shell=True)  
            for x in open(f'{name}.mlst').readlines():
                logging.info(f"MLST Output: {x}")
                if int(x.split()[2]) == int(records[x.split()[0].replace('.fasta', '' )]['ST']):
                    logging.info(f"Matches ST for {name}")

    


# Call MLST 
def call_mlst(): 
    pass

# check versus original table. 

if __name__ == "__main__":
    main()