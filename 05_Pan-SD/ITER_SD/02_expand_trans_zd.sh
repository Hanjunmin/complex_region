#####将expand.out转移到浙大上面去跑
cd /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expan_zd/

fold="/home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expan_zd/"
while read line; do
  echo "$line"
  line=$(echo "$line" | tr -d '\r')
  mkdir ${fold}$line
  cd ${fold}$line
  cp /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/01_iter/${line}/expand.out ./
done < /home/jmhan/PANSDEND/DATA/allsam

zstd -r 02_expan_zd ##迭代每个文件夹去压缩
rm /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expan_zd/*/expand.out

rsync -av /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expan_zd  clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/draft/
cd /share/home/zhanglab/user/maoyafei/project/PANSDEND/END/
rsync -av  clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/draft/02_expan_zd  ./
zstd -d -r  02_expan_zd
/share/home/zhanglab/user/maoyafei/project/PANSDEND/Z.run.sh


###如果在当前服务器跑
sample_list="/home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/weno2" #"/home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expan_zd/samples_without_zst.txt"
cat "$sample_list" | xargs -I{} -P 5 bash -c '
    sample=$(echo "{}" | tr -d '\''\r\n'\'')
    cd "/home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/01_iter/${sample}" || exit
    snakemake -s /home/jmhan/pansdall/iter.smk --config iterdir="/home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/01_iter/${sample}/" --keep-going --keep-incomplete -j 25
    echo $sample
'
###把当前服务器上的传到浙大上
cd /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expand_zd_trans/
while read line; do
mkdir /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expand_zd_trans/$line
cd /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expand_zd_trans/$line
mkdir minimapspnow
cd minimapspnow
cp /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/01_iter/$line/minimapspnow/all1.end.statistics ./
done < /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expand_zd_trans/weno
rsync -av /home/jmhan/PANSDEND/02_SAMSD/02_02-ITERSD/02_expand_zd_trans  clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/draft/

cd /share/home/zhanglab/user/maoyafei/project/PANSDEND/END/02_expan_zd
rsync -av clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/draft/02_expand_zd_trans1/* ./
