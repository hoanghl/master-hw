# scontrol -o show node output fields
# scontrol -o show node | gawk -F ' |=' '{for (i=1;i<=NF;i+=2) printf "%s|",$i; printf "\n"}' 

#    1        2           3             4         5         6
# NodeName | Arch | CoresPerSocket | CPUAlloc | CPUErr | CPUTot |
#
#     7        8         9     10           11          12
# CPULoad | Features | Gres | NodeAddr | NodeHostName | OS |
#
#    13          14        15      16          17            18
# RealMemory | Sockets | Boards | State | ThreadsPerCore | TmpDisk |
#
#   19       20            21               22             23
# Weight | BootTime | SlurmdStartTime | CurrentWatts | LowestJoules |
#
#       24
# ConsumedJoules |


BEGIN {
    FS="[ =]+";
    fmax=24;
    for (i=1;i<=fmax;i++) f[i]=0;
    f[1]=f[3]=f[4]=f[6]=f[7]=f[14]=f[16]=f[17]=1;
    fwidth=10;
    form=sprintf("%%-%ds ",fwidth);
}

NR==1 {
    for (i=1;i<=2*fmax;i+=2) {
	j=(i+1)/2;
	h[j]=$i;
    }
    for (i=1;i<=fmax;i++) {
	if (f[i]>0) printf form,substr(h[i],1,fwidth);
    }
    printf "\n";
}

{
    for (i=1;i<=fmax;i++) {
	if (f[i]>0) printf form,substr($(2*i),1,fwidth);
    }
    printf "\n";
}



