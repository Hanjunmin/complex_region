import os
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import argparse

def find_a_txt_files():
    return glob.glob("*.a.txt")

def get_sample_name(a_txt_file):
    """从a.txt文件名提取样本名"""
    return os.path.basename(a_txt_file).replace('.a.txt', '')

def filter_pairs(sample):
    """执行filter_pairs步骤"""
    a_txt = f"{sample}.a.txt"
    b_txt = f"{sample}.b.txt"
    output_file = f"filt/{sample}.end.txt"
    
    print(f"Processing file: {a_txt} and {b_txt}")
    
    # 创建输出目录
    os.makedirs("filt", exist_ok=True)
    
    # 执行R脚本
    cmd = ["Rscript", scr, a_txt, b_txt, output_file]
    subprocess.run(cmd, check=True)
    
    return output_file

def bash_pairs(sample):
    """执行bash_pairs步骤"""
    input_file = f"filt/{sample}.end.txt"
    output_file = f"filt/{sample}.filt.pair"
    
    cmd = f"""
    cat {input_file} | tr ':' '\\t' | awk -F '\\t' 'OFS="\\t"{{gsub(/-/, "\\t", $2);gsub(/-/, "\\t", $4);gsub(/-/,"\\t",$6);print $0;}}' | awk '((!($1 == $4 && $5 > $2 && $5 < $3) && $10 >=0.05) || ($1 == $4 && $5 > $2 && $5 < $3 && $10 > 0.5)){{print $0}}' > {output_file}
    """
    
    subprocess.run(cmd, shell=True, check=True)
    return output_file

def process_individual_pairs(sample):
    """执行process_individual_pairs步骤"""
    input_file = f"filt/{sample}.filt.pair"
    output_file = f"filt/{sample}.processed.pair"
    
    cmd = f"cat {input_file} | awk 'OFS=\"\\t\"{{print $1,$2,$3,$4,$5,$6}}' > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    return output_file

def process_sample(sample):
    """处理单个样本的完整流程"""
    try:
        print(f"开始处理样本: {sample}")
        
        # 步骤1: filter_pairs
        filter_pairs(sample)
        
        # 步骤2: bash_pairs  
        bash_pairs(sample)
        
        # 步骤3: process_individual_pairs
        process_individual_pairs(sample)
        
        print(f"完成处理样本: {sample}")
        return f"filt/{sample}.processed.pair"
        
    except Exception as e:
        print(f"处理样本 {sample} 时出错: {e}")
        return None

def combine_processed_pairs(processed_files):
    """合并所有处理后的文件"""
    output_file = "filt/beforegroup.txt"
    
    with open(output_file, 'w') as outfile:
        for file in processed_files:
            if file and os.path.exists(file):
                with open(file, 'r') as infile:
                    outfile.write(infile.read())
    
    print(f"所有文件已合并到: {output_file}")
    return output_file

def main():
    parser = argparse.ArgumentParser(description="并行处理*a.txt文件")
    parser.add_argument("--threads", "-t", type=int, default=4, 
                       help="并行线程数 (默认: 4)")
    parser.add_argument("--script", "-s",
                       help="script for filter pairs")
    args = parser.parse_args()
    global scr
    scr=args.script
    # 创建输出目录
    os.makedirs("filt", exist_ok=True)
    
    # 寻找所有*a.txt文件
    a_txt_files = find_a_txt_files()
    
    if not a_txt_files:
        print("未找到任何*a.txt文件")
        return
    
    samples = [get_sample_name(f) for f in a_txt_files]
    print(f"找到 {len(samples)} 个样本: {samples}")
    
    # 并行处理所有样本
    with ProcessPoolExecutor(max_workers=20) as executor:
        processed_files = list(executor.map(process_sample, samples))
    
    # 合并结果
    combine_processed_pairs(processed_files)
    
    print("所有处理完成!")

if __name__ == "__main__":
    main()
