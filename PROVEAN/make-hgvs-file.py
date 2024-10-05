#!/usr/bin/env python3


seq = 'MNINDLIREIKNKDYTVKLSGTDSNSITQLIIRVNNDGNEYVISESENESIVEKFISAFKNGWNQEYEDEEEFYNDMQTITLKSELN'

for del_start in range(len(seq)):
    for del_end in range(del_start,len(seq)):
        if del_start == del_end:
            hgvs = f"{seq[del_start]}{del_start+1}del"
        else:
            hgvs = f"{seq[del_start]}{del_start+1}_{seq[del_end]}{del_end+1}del"
        print(hgvs)
        
        
