import os.path
from FragMGen import FragMGen
hap_len = 1000
ploidy = 2
reads_len = 2
alleles_num = 2
insert1 = 150
insert2 = 250
gen_num = 10
prefix = 'Data22/P'+str(ploidy)+'-A'+str(alleles_num)+'/'
for cov_rate in [30,40,50]:#
    for err_prob in [0.1,0.2,0.3]:#
        output_dir = 'sim-cov'+str(cov_rate)+'-err'+str(err_prob)+'/'
        if not os.path.isdir(prefix+output_dir):
            os.makedirs(prefix+output_dir)
        for i in range(gen_num):
            output_file = prefix+output_dir+'sim'+str(i)+'.txt'
            output_file_v = prefix+output_dir+'sim'+str(i)+'-phased'+'.txt'
            FragMGen(cov_rate, hap_len, ploidy, err_prob, reads_len, alleles_num, insert1, insert2, output_file, output_file_v)