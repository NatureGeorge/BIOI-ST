---
author: Zefeng Zhu
date: 2020-09-08 09:11:02 +0800
---

# 1-2 Ab initio prediction and structural modeling of whole genome genes

## Design

![fig](https://user-images.githubusercontent.com/43134199/92424345-adc60880-f1b6-11ea-958d-b9d2f2d699c3.png)

### Organism

1. 草履虫 (原生生物)
2. 秀丽线虫 (线虫)
3. 酵母菌 (真菌)
4. 拟南芥 (被子植物)

### Software

1. Augustus
2. Genscan
3. GeneMark-ES/ET

### Procedures (Demo)

### 1. 基因组数据准备

> 下载基因组序列(FASTA 格式). 以及相应的 GFF 格式注释文件

```bash
workDir="/home/student/s24/zeFengZhu/Gen/lab4"
fastaFile="$workDir/GCA_000977265.3_Sc_YJM1342_v1_genomic.fna"
gffFile="$workDir/GCA_000977265.3_Sc_YJM1342_v1_genomic.gff"
```

### 2. 从头基因预测软件的安装与测试

#### `Augustus`

可用在Anaconda环境下利用Bioconda安装augustus,但仅针对Linux-64/Mac OSX-64系统而不支持Windows系统。

> Bioconda is a distribution of bioinformatics software realized as a channel for the versatile Conda package manager. (https://anaconda.org/bioconda/augustus)

```bash
conda install -c bioconda augustus
```

#### `Genscan`

pass

#### `GeneMark-ES/ET`

pass


### 3. 全基因组的从头基因预测

#### 3.1 对上述每个能够进行从头预测基因的软件进行实验

选择Augustus进行实验。

```bash
augustus --gff3=on --outfile=Sc_augustus_out.gff3 --species=saccharomyces_cerevisiae_S288C $fastaFile
augustus --species=saccharomyces_cerevisiae_S288C --UTR=off --strand=both --sample=100 --keep_viterbi=true --alternatives-from-sampling=false --genemodel=partial /data/www/augpred/webdata/pred9dTmEkZ9/genome.fa --codingseq=on --exonnames=on

```

> 对参数存疑

#### 3.2 使用该软件对第 1 步准备的基因组序列进行基因预测分析，

保存 GFF 格式的预测结果，以及相应的多肽或 CDS 序列(FASTA 格式)
得到文件：

```bash
augustus.aa
augustus.cdsexons
augustus.codingseq
augustus.gbrowse
augustus.gff
augustus.gtf
```

### 4 从头基因预测结果的鉴别

> 注： 4.1步骤在实验三中已经完成，本次实验直接采用实验三文件，下述再一次记录相关步骤。

#### 4.1 已知蛋白序列

根据基因组序列的物种来源，从 UniProt 数据库搜索. 下载近缘物种所有已知蛋白序列(reviewed)

进入UniProt进行搜索：

```bash
# 搜索内容
taxonomy:fungi NOT "saccharomyces cerevisiae" AND reviewed:yes
```

下载fasta格式文件：

```bash
# 得到文件
uniprot-taxonomy_fungi+NOT+_saccharomyces+cerevisiae_+AND+reviewed_yes.fasta
# 重命名为
protein.fasta
# 设置路径
unpFastaFile="$workDir/protein.fasta"
```

#### 4.2 创建本地 BLAST 数据库

使用 makeblastdb 程序，对上述 FASTA 格式的蛋白质序列进行处理，建立本地 BLAST 数据库

```bash
makeblastdb -in $unpFastaFile -input_type fasta -title uniprot_protein -dbtype prot -out uniprot_protein
```

输出如下：

```bash
Building a new DB, current time: 11/14/2019 17:06:57
New DB name:   /home/student/s24/zeFengZhu/Gen/lab4/uniprot_protein
New DB title:  uniprot_protein
Sequence type: Protein
Keep Linkouts: T
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 24905 sequences in 2.86177 seconds.
```

#### 4.3 从GFF文档中提取FASTA序列

##### GFF中序列格式范例

```bash
# start gene g1
...
# protein sequence = [MVKLTSIAAGQITSSITSSRPIITPFYPSNGTSVISSSVISSSVISSSVTSSL...
# SIFSESS...
...
# TTEITKQTTETTKQTTETTKQTTVVTIFSCESDVCSKTASPAIVSTSTATINDVTTEYTTWCPISTTESRQQT...]
# end gene g1
```

可以看到，FASTA序列记录的模式可以总结为:

```python
startwith = "# protein sequence = \[([A-z]+)" # # coding sequence = \[([a-z]+)
content = "([A-z]+)"
endwith = "# ([A-z]+)]"
endkey = "end gene ([A-z0-9]+)"
```

##### 提取序列函数

```python
def ExtractSeqFromGFF3(text, startwith=r"# protein sequence = \[([A-z]+)", content="([A-z]+)", endwith="([A-z]+)]", endKey="end gene ([A-z0-9]+)"):
    assert isinstance(text, (Iterable, Iterator)), "Invalid Object"

    startwith, content, endwith, endKey = (re.compile(i) for i in (startwith, content, endwith, endKey))
    flag, endToken, seq = 0, 0, ""

    for line in text:
        startToken = startwith.search(line)
        if startToken is not None:
            flag = 1
            seq += startToken.group(1)
        if flag:
            endToken = endwith.search(line)
            if endToken is not None:
                flag = 0
            if startToken is None:
                seq += content.search(line).group(1)
        elif endToken is not None:
            key = endKey.search(line)
            if key is not None:
                yield key.group(1), seq
                endToken, seq = 0, ""


def toFASTA(name, seq):
    return ">{name}\n{seq}\n".format(name=name, seq=seq)


def script(inPath, outPath, mode):
    with open(inPath, "rt") as inFile:
        with open(outPath, "wt") as outFile:
            if mode == "gene":
                g = ExtractSeqFromGFF3(inFile, startwith=r"# coding sequence = \[([a-z]+)")
            else:
                g = ExtractSeqFromGFF3(inFile)
            for name, seq in g:
                outFile.write(toFASTA(name, seq[:-1]))


if __name__ == "__main__":
    script("Sc_augustus_out.gff3", "augustus_gene.fasta", "gene")
    script("Sc_augustus_out.gff3", "augustus_protein.fasta", "protein")
```

提取出预测基因序列文件：```augustus_gene.fasta```; 蛋白序列文件：```augustus_protein.fasta```

#### 4.4 使用合适的 blast 程序对该预测基因与已知蛋白序列进行比对,以此来鉴别从头预测出来的基因

> 只保留打分最高的一条结果，由```max_target_seqs```指定

```bash
nohup blastx -query ./augustus_gene.fasta -db uniprot_protein -out ./Sc_blastx_gene_results.outfmt6 -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 10 > nohup_blastx_gene.out &
# 或
nohup blastp -query ./augustus_protein.fasta -db uniprot_protein -out ./Sc_blastp_gene_results.outfmt6 -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 10 > nohup_blastp_gene.out &
```


#### 4.5 把 4.4 结果合并到 3.2 获得的 GFF 格式结果中

使用相似性蛋白的缩写名称，替换原来预测的基因名称(模仿实验 3 格式转换后的 GFF 文档)，保存为一个新的 GFF 文档。

```bash
perl blast92gff3.pl Sc_blastx_gene_results.outfmt6 > Sc_blastx_gene_results.gff
# Summary of HSPs saved
# ALL saved = 4498
# other saved = 4498
```

##### 信息整合脚本

```py
def getMapping(filePath):
    dfrm_gff = pd.read_csv(filePath, sep="\t", header=None, skiprows=1)
    unp_pattern = re.compile("sp:([A-z0-9_\|]+)")
    gene_pattern = re.compile("Target=([A-z0-9]+)")

    di = {}

    for index in dfrm_gff.index:
        gene = gene_pattern.search(dfrm_gff.loc[index, 8]).group(1)
        di[gene] = unp_pattern.search(dfrm_gff.loc[index, 0]).group(1)

    return di

def updateGFF(inPath, outPath, di):
    with open(inPath, "rt") as inFile:
        with open(outPath, "wt") as outFile:
            startwith = re.compile("# start gene ([A-z0-9]+)")
            flag = 0
            for line in inFile:
                startToken = startwith.search(line)
                if startToken is not None:
                    flag = 1
                    key = startToken.group(1)
                    outFile.write(line)
                    continue
                if flag:
                    line = line[:-1] + ";%s\n" % di.get(key, "")
                    flag = 0
                outFile.write(line)


di = getMapping("Sc_blastx_gene_results.gff")
addNota("Sc_augustus_out.gff3", "augustus_addNota.gff", di)
```

变量```di```即可具体查看匹配上基因的情况，进行统计。

得到加入了蛋白缩写名称的augustus结果```augustus_addNota.gff```。

### 5. 从头预测结果的评估

#### 5.1 gffcompare对比

使用 gffcompare 工具把第 4 步结果与 1.1 步原始 GFF 数据以及实验 3 结果进行比较， 查看结果，并分析它们之间的异同之处。

```bash
gffcompare -V -r $gffFile ./augustus_addNota.gff -o ./Sc_augustus_out_addNota
```

得到下列文件：

```bash
Sc_augustus_out_addNota.annotated.gtf
Sc_augustus_out_addNota.augustus_addNota.gff.refmap
Sc_augustus_out_addNota.augustus_addNota.gff.tmap
Sc_augustus_out_addNota.loci
Sc_augustus_out_addNota.stats
Sc_augustus_out_addNota.tracking
```

#### 5.2 gffcompare结果解析

提取```Sc_augustus_out_addNota.stats```内容如下：

|| Sensitivity | Precision  |
-|-|-
Base level|    99.5     |    95.6    |
Exon level|    93.7     |    88.8    |
Intron level|    72.4     |    50.1    |
Intron chain level|    71.9     |    52.4    |
Transcript level|    95.4     |    92.0    |
Locus level|    95.5     |    92.0    |
-|-|-
Matching intron chains|164|
Matching transcripts|4956|
Matching loci|4955|
-|-|-
Missed exons|48/5434|(0.9%)
Novel exons|333/5730|(5.8%)
Missed exons|21/239|(8.8%)
Missed exons|125/345|(36.2%)
Missed exons|0/5187|(0.0%)
Missed exons|251/5385|(4.7%)


![fig](https://github.com/NatureGeorge/ProBioinformatics/raw/master/Omics/Genomics/final/figs/GFF.png "Fig of GFF")

可以知道，GFF文档相当于在检验BLAST比对找到的UniProt与基因组GFF注释文档里的UniProt的结果接近程度。且视基因组GFF注释文档内容皆为真。

* Base Level：在相同坐标上报告的外显子碱基的数目情况
  * Sensitivity：高达99.5%，说明BLAST结果在该水平上结果找到了绝大部分基因组注释文档中的内容，极少数碱基根本没有被任何预测的转录本(transfrags)外显子所覆盖
  * Precision：达95.6%，说明BLAST结果在该水平上有一小部分(4.4%)碱基被预测的转录本外显子覆盖但未被任何参考转录本外显子覆盖
* Exon level：两文件基因组上的外显子间隔交集情况
  * 可以看到预测基因结果的外显子与基因组注释文档的外显子边界有一定小差异
* Intron level：内含子间隔
  * 预测基因的内含子边界有不少与基因组注释文档存在差异，且错误预测了更多内含子，Precision仅50.1%
* Transcript level：预测转录本与参考转录本间的匹配情况
  * 转录水平是匹配良好，但也有少数"误差"
* Locus level：观察到的基因座(外显子重叠的转录物簇)与构建的参考基因座的相似匹配情况
  * 基因座位置也匹配良好，但也有少数"误差"

## Division of tasks

pass