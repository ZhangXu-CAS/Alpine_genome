# -*- coding: utf-8 -*-
# callCodeML 
# Verson: 2.0
# Date: 2021.11.15
# This is used to call Codeml to calculate the positive selection of gene in Branch-Sites Model.(Require codon sequence)
# Usage: python3 callCodeml.py seqDir treeFile
# Refactoring code and running logic to support multiprocess

import re
import sys
import os
from scipy.stats import chi2
import time


def help_info():
    print("This script is used to call codeML.\nUsage: python3 callCodeml.py Dir treeFile\n")


def getArgs():
    global workDir 
    global tree_file
    global dir_path 
    
    start_time = time.strftime('%Y%m%d-%H%M%S') # to keep the working dircitory name
    workDir = f"{os.getcwd()}/WorkingDir_{start_time}"

    args = sys.argv[1:]
    if len(args) == 0:
        help_info()
    if os.path.isdir(args[0]) is True:
        dir_path = os.path.abspath(args[0])  # path of seqences
        seqs = os.listdir(dir_path)
        seqs = [x for x in seqs if not str(x).startswith('.') ] # remove hided files
        tree_file = os.path.abspath(args[1])
        print(f"The drictory is: {dir_path}\nThe treeFile is: {tree_file}")
        #os.chdir(os.path.abspath(args[0]))
        return (dir_path, seqs, tree_file)
    else:
        help_info()

def create_dir(name):
    try:
        os.mkdir(f'{workDir}/')
    except:pass
    try:
        os.mkdir(f'{workDir}/{name}')
    except:pass
    try:
        os.mkdir(f'{workDir}/{name}/null')
        os.mkdir(f'{workDir}/{name}/alte')
    except:pass


def creat_ctl(name, tree_file):


       
    ctl_null = f'''
    seqfile = {dir_path}/{name}
    treefile = {tree_file}
    outfile = {workDir}/{name}/null/{name}_null.res

    noisy = 9   
    verbose = 0   
    runmode = 0  

    seqtype = 1   
    CodonFreq = 2   
    clock = 0   
    aaDist = 0
    model = 2

    NSsites = 2   
    icode = 0   
    Mgene = 0

    fix_kappa = 0   
    kappa = 2   
    fix_omega = 1   
    omega = 1   

    fix_alpha = 1   
    alpha = .0  
    Malpha = 0   
    ncatG = 3   

    getSE = 0   
    RateAncestor = 0   

    fix_blength = 0  
    method = 0   
    Small_Diff = .45e-6
    cleandata = 1
    '''
    ctl_alte = ctl_null.replace(f'{workDir}/{name}/null/{name}_null.res',  f'{workDir}/{name}/alte/{name}_alte.res')\
        .replace('fix_omega = 1', 'fix_omega = 0')\
        .replace('omega = 1', 'omega = 2')

    with open(f'{workDir}/{name}/null/null_profile.ctl', 'w') as f:
        f.write(ctl_null)
    
    with open(f'{workDir}/{name}/alte/alte_profile.ctl', 'w') as f:
        f.write(ctl_alte)


def run_create(dir_path, seqs, tree_file):
    
    for i in (seqs):
        
        name= i.split('/')[-1]
        create_dir(name)
        creat_ctl(name, tree_file)
    
    with open(f'{workDir}/codeml.sh', 'w') as f:
        f.write('''
#!/bin/bash

list=`ls -1 |grep OG`
count=`ls -1 |grep OG|wc|awk '{print $1}'`
echo "Total number: $count"
dir=`pwd`
echo "Start Now..."

ST=$(date +%s)
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1 # 创建有名管道
exec 5<>/tmp/fd1 # 创建文件描述符，以可读（<）可写（>）的方式关联管道文件，文件描述符5拥有有名管道文件的所有特性
rm -rf /tmp/fd1 # 文件描述符关联后拥有管道的所有特性，所有删除管道
NUM=11 # 获取输入的并发数

for (( i=1;i<=${NUM};i++ ))
do
  echo >&5 # &5表示引用文件描述符5，往里面放置一个令牌
done

for i in $list
do
  read -u5
  {
    cd $dir/$i/null
    codeml ./null_profile.ctl  >log.txt
    cd $dir/$i/alte
    codeml ./alte_profile.ctl  >log.txt
    sleep 1 # 模拟程序、命令
    echo >&5 # 执行完把令牌放回管道
  }& # 把循环体放入后台运行，相当于是起一个独立的线程，在此处的作用就是相当于起来10个并发
done
wait # wait命令的意思是，等待（wait命令）上面的命令（放入后台的）都执行完毕了再往下执行，通常可以和&搭配使用

ET=$(date +%s)
TIME=$(( ${ET} - ${ST} ))
echo "time: ${TIME}"
''')
        
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} Files created.\n")
    


    
def run_codeml():
    os.chdir(workDir)
    os.system('bash codeml.sh')
    print(f"\n{time.strftime('%Y-%m-%d %H:%M:%S')} Codeml Done...\n")


def run_stat(name):

    def get_res(path):
        with open(path, "r") as f:
            t =  f.read()
            kappa = re.findall(r"kappa.+", t)[0].split(' ')[-1].replace('\n', '')
            ls = re.findall(r'lnL.+', t)[0].split(' ')
            ls=[x for x in ls if x!='']
            lnL = ls[4]
            np = ls[3].replace("):","")
            return([lnL, np, kappa])

    path1 = f'{workDir}/{name}/alte/{name}_alte.res'
    path2 = f'{workDir}/{name}/null/{name}_null.res'
    
    try:
        res_alte = get_res(path1)
        res_null = get_res(path2)
    
        lnL0 = float(res_null[0])
        lnL1 = float(res_alte[0])
        np0 = int(res_null[1])
        np1 = int(res_alte[1])
        kappa = float(res_alte[2])

        lnl = abs(lnL0 - lnL1 )*2
        np = abs(np0 - np1)
        pvalue = 1 - chi2.cdf(lnl, np)


        #if pvalue <= 0.05:
        with open(f'{workDir}/result.tsv', 'a') as f:
            f.write(f"{name}\t{lnL0}\t{lnL1}\t{np0}\t{np1}\t{kappa}\t{pvalue}\n")
    except:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {name} failed to stats.")


def main():   
    
    arges = getArgs()
    run_create(arges[0], arges[1], arges[2])
    run_codeml()

    for i in  os.listdir(workDir):
        run_stat(i)


 

if __name__ == "__main__":
    t0 = time.time()
    
    main()
   
    print(f'Total time used: {time.time() - t0} S\n')

