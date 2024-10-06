#!/usr/bin/env python3

seq = 'MNINDLIREIKNKDYTVKLSGTDSNSITQLIIRVNNDGNEYVISESENESIVEKFISAFKNGWNQEYEDEEEFYNDMQTITLKSELN'
print("wt_seq,mut_seq,start_pos")

for del_start in range(len(seq)):
    for del_size in range(1,11):
        del_end = del_start + del_size
        if del_end > len(seq):
            continue
        mut_seq = seq[:del_start] + seq[del_end:]
        print(f"{seq},{mut_seq},{del_start+1}")
        #print(del_start+1, del_size, mut_seq)


