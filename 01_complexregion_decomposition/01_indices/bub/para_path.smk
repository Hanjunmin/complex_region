with open("al.subregional.1M") as f:
    lines = [line.strip().split('\t') for line in f.readlines()]


FOLD="01_GFA/"
pairs = [(line[0], (line[1] + ":" + line[2] + "-" + line[3]),line[1]) for line in lines]  # example：[(A1, B1), (A2, B2), ...]

output_files = [FOLD+f"{C}/{A}/{B}.OK" for A, B,C in pairs]


rule all:
    input:
        output_files

rule generate_file:
    input:
       
    output:
        FOLD+"{C}/{A}/{B}.OK"
    shell:
        """
        subfold={wildcards.A}
        segm={wildcards.B}
        chr=$(echo $segm |tr ":" "\t" |tr "-" "\t"| cut -f1)
        start=$(echo $segm |tr ":" "\t" |tr "-" "\t" | cut -f2)
        end=$(echo $segm |tr ":" "\t" |tr "-" "\t"| cut -f3)
        len=$(echo -e "$end-$start" | bc)
        echo $len
        cd {FOLD}$chr/$subfold
        echo $(pwd)
        cat ${{segm}}.gfa |grep '^S' |cut -f 2-3 >"temp."${{segm}}.S.start
        cat ${{segm}}.path.dou.filt  >"temp."${{segm}}.allpath
        python POLY_complex_dupseg_quick.py --sim 0.9 --mb 50 --dell 200  --region ${{segm}} 2>${{segm}}.logerr.txt

        cat ${{segm}}.gfa >end.gfa
        cat ${{segm}}.gfa NODE.P.${{segm}}data >${{segm}}.test.gfa  || true 
        cat <(cat ${{segm}}.gfa |head -n 1) S.${{segm}}.data L.${{segm}}.data P.${{segm}}.data >${{segm}}.sim.gfa  || true
        python /home/jmhan/COMPLEXFEATURE/09_bubble_only/01_GFA/simnode_inte.py --reg ${{segm}}
        python /home/jmhan/COMPLEXFEATURE/09_bubble_only/01_GFA/poly_ratio.py --reg ${{segm}}
        [ -f "temp."${{segm}}.S.start ] && rm "temp."${{segm}}.S.start
        [ -f "temp."${{segm}}.allpath ] && rm "temp."${{segm}}.allpath
        if [ -f polycomplex${{segm}}.txt ]; then
            touch {output}
        fi
        """