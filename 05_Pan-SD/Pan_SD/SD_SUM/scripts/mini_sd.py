import os
import subprocess
from pathlib import Path

# ---------------------------------------
# 工具函数
# ---------------------------------------

def run_cmd(cmd, shell=True):
    print(f"[CMD] {cmd}")
    subprocess.run(["bash", "-c", cmd], check=True)
    

def file_exists(path):
    return Path(path).is_file() and os.path.getsize(path) > 0

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

# ---------------------------------------
# 1. 前期准备
# ---------------------------------------

def prepare_dirs(base_dir, line):
    line_dir = Path(f"{base_dir}/02_refSD/sam_cor/2_ALIGN/{line}")
    t2t_dir = line_dir / "t2tsdminimap"
    if t2t_dir.exists():
        run_cmd(f"rm -r {t2t_dir}")
    ensure_dir(t2t_dir)
    return line_dir, t2t_dir


# ---------------------------------------
# 2. 生成中间文件 (.sd/.can → tsv)
# ---------------------------------------

def generate_sd_can_tables(line):
    cmds = [
        f"less -S {line}.sd | tail -n +2 | awk 'OFS=\"\\t\"{{print $7\"#\"$9,$10+$3,$10+$4}}' > allt2tsdcoorend.tsv",
        f"less -S {line}.sd | tail -n +2 | awk 'OFS=\"\\t\"{{print $9,$10+$3,$10+$4}}' > allt2tsdcoorend_true.tsv",
        f"less -S {line}.can | tail -n +2 | awk 'OFS=\"\\t\"{{print $7\"#\"$9,$10+$3,$10+$4}}' > candidateall.tsv",
        f"less -S {line}.can | tail -n +2 | awk 'OFS=\"\\t\"{{print $9,$10+$3,$10+$4}}' > candidateall_true.tsv"
    ]
    for cmd in cmds:
        run_cmd(cmd)


# ---------------------------------------
# 3. 调用 Rscript T2Tmatch
# ---------------------------------------

def run_t2tmatch(script, base_dir, line):
    cmd = f"Rscript {script}ref_SD/sam_sd/scripts/T2Tmatch.r {base_dir}/02_refSD/sam_cor/1_SD_CANPATH/sd.edge {line}.sd"
    run_cmd(cmd)


# ---------------------------------------
# 4. 生成 reprocess.pair、reprocess.pairall 等文件
# ---------------------------------------

def generate_reprocess_files(line):
    cmds = [
        r"""bash -c 'cat <(cat process.txt |grep -v "nosam" |tail -n +2 |cut -f1,2 |tr "\t" "\n" |sort |uniq) <(cat end.txt |tail -n +2|cut -f1,2 |tr "\t" "\n") > reprocess.pair'""",

        r"""bash -c 'cat <(cat process.txt|grep -v "nosam" |tail -n +2 |cut -f1,2) <(cat end.txt |tail -n +2 |cut -f1,2) > reprocess.pairall'""",

        r"""bash -c 'cat process.txt |grep "nosam" > processnosam.txt'""",
        f"""bash -c 'grep -Ff <(cut -f1 reprocess.pair) {line}.sd > {line}.multicoordi.inte.tsv'"""
    ]
    for cmd in cmds:
        subprocess.run(cmd, shell=True, check=False)


# ---------------------------------------
# 5. 主条件分支: 是否存在 multicoordi 文件
# ---------------------------------------

def handle_multi_coordinate(script, line):
    multi = f"{line}.multicoordi.inte.tsv"
    if file_exists(multi):
        print("[INFO] 多坐标情况：处理 pair 逻辑")
        cmds = [
            f'cat {multi} | awk \'OFS="\\t"{{print $2,$7"@"$9"@"$10"@"$11"@"$3"@"$4}}\' > cor.id',
            f"Rscript {script}ref_SD/sam_sd/scripts/intepair.r",
            'less -S pairall |cut -f 3,4 |tr "@" "\\t" | awk \'OFS="\\t"{print $2,$3+$5,$3+$6,$8,$9+$11,$9+$12}\' | tail -n +2 > pairidall',
            'paste pairidall <(cat pairall | tail -n +2|cut -f 2) <(cat pairall|tail -n +2 |cut -f 1) > T2T2SAM',
            'less -S T2T2SAM |awk \'OFS="\\t"{print $1":"$2"-"$3,$7,$4":"$5"-"$6,$8}\' | awk \'OFS="\\t"{print $1,$2"\\n"$3,$4}\' | sort |uniq > T2T2SAM.pair'
        ]
        for cmd in cmds:
            run_cmd(cmd)
        return True
    else:
        print("[INFO] 无坐标情况：生成空文件")
        for f in ["T2T2SAM", "T2T2SAM.pair", "pairidall"]:
            Path(f).touch()
        for f in ["T2Thave.temp", "T2Tpolymor1.temp", "other_T2T", "all.end.statistics", "all.end.statistics.filt","T2Tpolymor1other.bed"]:
            (t2t_dir / f).touch()
        
        return False


# ---------------------------------------
# 6. minimap2 对齐
# ---------------------------------------

def run_minimap2(script, genomefa, line, t2t_dir):
    fa = list(Path(genomefa).glob(f"{line}*.fa"))[0]
    fa = str(fa)
    run_cmd(f"cat {fa} > {t2t_dir}/{line}.fa")
    run_cmd(f"samtools faidx {t2t_dir}/{line}.fa")
    all_paf = t2t_dir / "all.paf"
    if all_paf.exists():
        all_paf.unlink()
    pairidall = Path("../pairidall")
    if not file_exists(pairidall):
        print("[WARN] pairidall 不存在，跳过 minimap2")
        return
    with open(pairidall) as f:
        for line_pair in f:
            col1, col2, col3, col4, col5, col6 = line_pair.strip().split("\t")
            cmd = (
                f'minimap2 -c --eqx <(samtools faidx {fa} "{col1}:{col2}-{col3}") <(samtools faidx {fa} "{col4}:{col5}-{col6}") -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >> {all_paf}'
            )
            run_cmd(cmd)
    cmds = [f"python {script}ref_SD/sam_sd/scripts/paf2bed.py --paf {all_paf} --o all.end.statistics",
            "awk '$10>=0.9 && $11>=1000 && $12>=0.5' all.end.statistics > all.end.statistics.filt",
            "awk '$15>=0.8 && $16>=0.8' all.end.statistics.filt > T2Thave.temp",
            "awk '!($15>=0.8 && $16>=0.8)' all.end.statistics.filt > T2Tpolymor.bed",
            "awk '$15>=0.95 || $16>=0.95' T2Tpolymor.bed > T2Tpolymor1.temp",
            "awk '!($15>=0.95 || $16>=0.95)' T2Tpolymor.bed > T2Tpolymor1other.bed",
            "cut -f14 T2Tpolymor1other.bed > all_values.txt",
            f"""Rscript -e 'library(data.table); \
    data<-fread("all_values.txt"); \
    data2<-fread("../T2T2SAM.pair", header=FALSE); \
    colnames(data2)<-c("quesou","T2T"); \
    df <- merge(data, data2, by="quesou", all.x=TRUE); \
    fwrite(df,"other_T2T", sep="\\t")'"""
        ]
    for c in cmds:
        run_cmd(c)

# ---------------------------------------
# 7. Rscript 过滤、bedtools 等操作
# ---------------------------------------

#########################################
def handle_T2Thave(t2t_dir, script_dir):
    f = "T2Thave.temp"
    if file_exists(f):
        run(f"Rscript {script_dir}ref_SD/sam_sd/scripts/filter_sam_SD.r T2Thave.temp T2Thave.bed", cwd=t2t_dir)
    else:
        Path("T2Thave.bed").touch()
    f = "T2Tpolymor1.temp"
    if file_exists(f):
        run(f"Rscript {script_dir}ref_SD/sam_sd/scripts/filter_sam_SD.r T2Tpolymor1.temp T2Tpolymor1.bed", cwd=t2t_dir)
    else:
        Path("T2Tpolymor1.bed").touch()



# Part 3
#########################################
def handle_other_T2T(t2t_dir, line, script, base_rm):
    other = t2t_dir / "other_T2T"
    cmds = []
    if file_exists(other):
        cmds.append(
            f'less {other} | awk \'OFS="\\t"{{print $2,$1}}\' | '
            f'sort | uniq | tr ":" "\\t" | '
            f'awk \'BEGIN {{OFS="\\t"}} {{gsub("-", "\\t", $2); gsub("-", "\\t", $4); print}}\' '
            f'| grep -v "T2T" > {t2t_dir}/t2trm'
        )
        cmds.append(
            f'bedtools intersect -a {t2t_dir}/t2trm -b <(cut -f1-3 {base_rm}) -wa -wb > {t2t_dir}/groupin'
        )
        cmds.append(
            f'bedtools subtract -a {t2t_dir}/t2trm -b <(cut -f1-3 {base_rm}) -A | cut -f4-6 > {t2t_dir}/groupin.ana'
        )
        cmds.append(
            f'Rscript {script}ref_SD/sam_sd/scripts/expand_rmtrf.r {t2t_dir}/groupin {line}.fa.fai {t2t_dir}/rm_in 2'
        )
        cmds.append(
            f'cat {t2t_dir}/rm_in {t2t_dir}/groupin.ana > {t2t_dir}/rm_in.end'
        )
        cmds.append(
            f'bedtools coverage -a {t2t_dir}/rm_in.end -b ../RMTRF/RMTRFSAT.region | '
            f'awk \'$7!=1\' > {t2t_dir}/rmin.in'
        )
        cmds.append(
            f'bash {script}ref_SD/sam_sd/scripts/run_rm_trf.sh {t2t_dir}/rmin.in {line} /home/Public/software/RepeatMasker-4.1.4/ {script}ref_SD/sam_sd/scripts/'
        )
        for c in cmds:
            run_cmd(c)
    else: 
        for f in ["rm.out", "rm_trf.s.bed", "rmin.in"]:
            (t2t_dir/f).touch()

    


#########################################
# Part 4
#########################################
def handle_rm_trf(t2t_dir, script_dir, line, base_dir):
    cmds = [
        f"bash -c 'cat <(cat ../RMTRF/RMSAT.sat | awk \"OFS=\\\"\\t\\\"{{print $0,\\\"Satellite\\\"}}\") "
        f"<(cat rm.out | cut -f 1-3,5) > allsat'",
        "cat rm_trf.s.bed ../RMTRF/RMTRF.rmtrf > allrmtrf",
        f"bash {script_dir}ref_SD/sam_sd/scripts/transe.sh allrmtrf allrmtrf.trane",
        f"bash {script_dir}ref_SD/sam_sd/scripts/endfilt.sh allsat allrmtrf.trane T2Tpolymor1other.bed filin.bedend"]
    for c in cmds:
        run_cmd(c)

def handle_rm_trf1(script_dir, line, base_dir):
    cmds = [
        f"cp {base_dir}/02_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/filin.bedend {base_dir}/02_refSD/sam_cor/2_ALIGN/{line}/filin.bedend",
        "cat process.txt end.txt | grep -v 'sam' > sim.txt || true",
        "pwd",
        f"cat {line}.sd | tail -n +2 | awk 'OFS=\"\\t\"{{print $2,$9\":\"$10+$3\"-\"$10+$4,$5}}' > pos.dic",
        f"Rscript {script_dir}SD_SUM/scripts/t2t_mer_al.r {line} {base_dir}/02_refSD/sam_cor/1_SD_CANPATH/sd.edge"
    ]
    for c in cmds:
        run_cmd(c)

# ---------------------------------------
# 8. 主入口
# ---------------------------------------
# base_dir="/share/home/zhanglab/user/maoyafei/complex_github/RUN_TEST/RUN/"
# line="C002-CHA-E02-Mat"
# script="/share/home/zhanglab/user/maoyafei/complex_github/05_Pan-SD/Pan_SD/"
# genomefa="/share/home/zhanglab/user/maoyafei/complex_github/05_Pan-SD/example/FA/"
# base_rm="/share/home/zhanglab/user/maoyafei/complex_github/RUN_TEST/RUN//01_intersim/CHM13v2.fa_rm.bed"

def main(base_dir, genomefa, script, line,base_rm):
    global t2t_dir
    line_dir, t2t_dir = prepare_dirs(base_dir, line)
    os.chdir(line_dir)
    
    generate_sd_can_tables(line)
    run_t2tmatch(script, base_dir, line)
    generate_reprocess_files(line)
    has_multi = handle_multi_coordinate(script, line)
    os.chdir(t2t_dir)
    run_minimap2(script, genomefa, line, t2t_dir)
    handle_T2Thave(t2t_dir, script)
    handle_other_T2T(t2t_dir, line, script, base_rm)
    handle_rm_trf(t2t_dir, script, line, base_dir)
    os.chdir(line_dir)
    handle_rm_trf1(script, line, base_dir)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--genomefa", required=True)
    parser.add_argument("--script", required=True)
    parser.add_argument("--line", required=True)
    parser.add_argument("--base_rm", required=True)

    args = parser.parse_args()

    main(args.base_dir, args.genomefa, args.script,  args.line, args.base_rm)
