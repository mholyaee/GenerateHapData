import numpy as np

def FragMGen(cov_rate, hap_len, ploidy, err_prob, reads_len, alleles_num, insert1, insert2, output_file, output_file_v):
    # Finding number of reads
    seq_cov = cov_rate*ploidy
    read_num = np.int(np.ceil((hap_len-reads_len + 1)/reads_len*seq_cov)+1)
    # constructing V - Haplotype sequence matrix
    V = np.random.randint(1, alleles_num+1, size  = [ploidy,hap_len]) 
    # removing homozygous sites
    tt=0
    if ploidy == 2:
        V[1,:] = -(V[0,:]-3)
    else:
        for ii in range(hap_len):
            if np.size(np.unique(V[:,ii])) == 1:
                tt = tt+1
                temp1=np.random.randint(ploidy)
                temp = np.mod(np.random.randint(alleles_num-1) + np.mod(V[temp1,ii],alleles_num),alleles_num)
                V[temp1,ii] = (temp + alleles_num*(temp == 0))
    # constructing U - Haplotype membership matrix
    U = np.random.randint(ploidy, size = [read_num])
    I = np.eye(ploidy)
    U = I[U]
    # constructing M - True SNP Fragment Matrix
    M = (np.dot(U,V)).astype(int)-1
    V = np.transpose(V)-1
    # constructing R - Partially observed-noisy SNP Fragment Matrix
    fpo=open(output_file,'w')
    fpo.write(str(read_num)+'\n')
    fpo.write(str(hap_len))
    for ii in range(read_num):
        # Masking structure - Paired-end reads
        r1 = np.random.poisson(lam = reads_len)
        r2 = np.random.poisson(lam = reads_len)
        # Check whether the read is informative
        if r1 + r2 < 2:
            continue
        else:
            fpo.write('\n')
            d = np.random.randint(1,insert1+1) + insert2
            st_point = np.random.randint(hap_len-r1-r2-d+1)
            # Effect of error for first read block 
            if r1==0:
                Mr1 = np.array([])
            else:
                Mr1 = M[ii,np.array(range(r1))+st_point]
                temp = np.mod(np.random.randint(1,alleles_num,size = [np.size(Mr1)]) + np.mod(Mr1,alleles_num),alleles_num)
                rn = (np.random.random((np.size(Mr1),)) >= err_prob)
                Mr1 = (rn*Mr1 + np.logical_not(rn)*(temp)).astype(int) 
            # Effect of error for second read block 
            if r2==0:
                Mr2 = np.array([])
            else:
                Mr2 = M[ii,np.array(range(r2))+st_point+r1+d]
                temp = np.mod(np.random.randint(1,alleles_num,size = [np.size(Mr2)]) + np.mod(Mr2,alleles_num),alleles_num)
                rn = (np.random.random((np.size(Mr2),)) >= err_prob)
                Mr2 = (rn*Mr2 + np.logical_not(rn)*(temp)).astype(int)
            if r1 == 0:
                fpo.write('1')
                fpo.write(' RID_'+str(ii)+' '+str(st_point+r1+d)+' '+''.join(str(e2) for e2 in Mr2))
            elif r2 == 0:
                fpo.write('1')
                fpo.write(' RID_'+str(ii)+' '+str(st_point)+' '+''.join(str(e1) for e1 in Mr1))
            else:
                fpo.write('2')
                fpo.write(' RID_'+str(ii)+' '+str(st_point)+' '+''.join(str(e1) for e1 in Mr1)+' '+str(st_point+r1+d)+' '+''.join(str(e2) for e2 in Mr2))
    np.savetxt(output_file_v,V,fmt='%d')