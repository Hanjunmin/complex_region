#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor
import sys

def load_config():
    """加载配置参数"""
    # 这里可以根据需要从配置文件、环境变量或命令行参数获取
    config = {
        'base_dir': None,
        'fa': None
    }
    return config

def find_working_folders(base_dir):
    """寻找需要处理的文件夹"""
    all_folders = glob.glob(os.path.join(base_dir, "*/"))
    print(f"找到所有文件夹: {all_folders}")
    
    # 排除不需要处理的文件夹
    exclude_folders = ['minimapsplit', 'minimapspnow', 'endfilt', 'split']
    folder_names = [
        os.path.basename(f.rstrip('/')) 
        for f in all_folders 
        if os.path.basename(f.rstrip('/')) not in exclude_folders
    ]
    print(f"需要处理的文件夹: {folder_names}")
    return folder_names


def process_single_folder(folder_name, base_dir, fa_file):
    try:
        txt_file = folder_name+"/"+folder_name+".txt"
        clu1_file =folder_name+"/"+folder_name+".clu1"
        log_file = os.path.join(base_dir, folder_name, f"{folder_name}.log")
        cmd = [
            "Rscript", 
            scr,
            txt_file,
            clu1_file,
            folder_name,
            base_dir
        ]
        
        print(f"执行命令: {' '.join(cmd)}")
        
        # 运行命令并重定向输出到日志文件
        with open(log_file, 'w') as log_f:
            result = subprocess.run(cmd, stdout=log_f, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            print(f"成功处理文件夹: {folder_name}")
            return True
        else:
            print(f"处理文件夹 {folder_name} 失败，错误: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"处理文件夹 {folder_name} 时发生异常: {e}")
        return False

def create_allclus_file(base_dir):
    """创建allclus文件（如果需要）"""
    allclus_file = os.path.join(base_dir, "allclus")
    try:
        with open(allclus_file, 'w') as f:
            f.write("# All clusters processed\n")
        print(f"创建allclus文件: {allclus_file}")
        return True
    except Exception as e:
        print(f"创建allclus文件失败: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="并行处理文件夹中的R脚本")
    parser.add_argument("--base_dir", "-b", required=True, 
                       help="基础目录路径")
    parser.add_argument("--fa", "-f", required=True,
                       help="FASTA文件路径")
    parser.add_argument("--threads", "-t", type=int, default=10,
                       help="并行线程数 (默认: 10)")
    parser.add_argument("--script", "-s",
                       help="script for filter pairs")
    
    args = parser.parse_args()
    global scr
    scr=args.script
    if not os.path.exists(args.base_dir):
        print(f"错误: 基础目录不存在: {args.base_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.fa):
        print(f"错误: FASTA文件不存在: {args.fa}")
        sys.exit(1)
    
    # 查找需要处理的文件夹
    folder_names = find_working_folders(args.base_dir)
    
    if not folder_names:
        print("没有找到需要处理的文件夹")
        sys.exit(0)
    
    print(f"找到 {len(folder_names)} 个需要处理的文件夹")
    
    
    # 并行处理所有文件夹
    successful_folders = []
    failed_folders = []
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        # 提交所有任务
        future_to_folder = {
            executor.submit(process_single_folder, folder, args.base_dir, args.fa): folder 
            for folder in folder_names
        }
        
        # 收集结果
        for future in future_to_folder:
            folder_name = future_to_folder[future]
            try:
                success = future.result()
                if success:
                    successful_folders.append(folder_name)
                else:
                    failed_folders.append(folder_name)
            except Exception as e:
                print(f"处理文件夹 {folder_name} 时发生未预期错误: {e}")
                failed_folders.append(folder_name)
    
    base_dir=args.base_dir
    cmd = f'find {base_dir} -name "*_largecluster.txt" -exec cat {{}} + > largecluster.bed || touch largecluster.bed; ' \
          f'find {base_dir} -name "*_cluster.bed" -exec cat {{}} + > cluster.txt || touch cluster.txt; ' \
          f'cat largecluster.bed cluster.txt > {base_dir}/allclus'
    subprocess.run(cmd, shell=True, check=True)


if __name__ == "__main__":
    main()