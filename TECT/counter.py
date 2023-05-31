import os
import re
import pysam
import copy
import pickle
import pandas as pd
import gzip

from TECT.tools import *
from datetime import datetime


class Counter:
    def __init__(self, index, out):  # 使用rmsk，gff进行初始化；包含out路径
        self.out_path = out  # 输出的路径

        if self.out_path[-1] != '/':
            self.out_path += '/'
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)

        self.index_path = index

        self.te_annotation = {}  # 每条染色体->{每个基因区域的起始位置->(end_pos,b_rc,sub_type,family,csn_start,csn_end)}
        self.interval_tree_te = {}  # 构建间隔树，字典key为chr，value里存树

        self.interval_tree_exon = {}

        self.TE_counter = {}  # 记录TE数量
        self.df = {}

        self.set_f_bc=set()

        # 载入index文件
        self.load_index_file(self.index_path)

    def get_filtered_barcodes(self,path):
        path=path.rstrip('/')
        filename=os.path.basename(path)
        suffix=filename.split('.')[-1]
        if suffix=='gz':
            f=gzip.open(path,'rt')
            for line in f:
                f_bc=line.rstrip()
                if '-' in f_bc:
                    f_bc=f_bc.split('-')[0]
                self.set_f_bc.add(f_bc)
            f.close()
        else:
            with open(path, 'r') as f:
                for line in f:
                    f_bc=line.rstrip()
                    if '-' in f_bc:
                        f_bc=f_bc.split('-')[0]
                    self.set_f_bc.add(f_bc)

    def load_index_file(self, index):
        print('[%s] Start the TECT_count module.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        print('')
        print('[%s] Loading the index file...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        with open(index, 'rb') as f_index:
            idx_file = pickle.load(f_index)
        self.te_annotation = idx_file['te']
        self.interval_tree_te = idx_file['te_tree']
        self.interval_tree_exon = idx_file['exon_tree']
        self.TE_counter = idx_file['counter']
        print('[%s] The index file has been load.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def new_df(self):
        subs = []
        for fam in self.TE_counter:
            for sub in self.TE_counter[fam]:
                subs.append(sub)
        subs.extend(['sum_exon', 'sum_map'])
        # print(subs)
        self.df = pd.DataFrame(columns=subs)

    def is_within_repeat_region_interval_tree(self, chrm, pos):
        if chrm not in self.interval_tree_te:
            return False, -1
        tmp_tree = self.interval_tree_te[chrm]
        set_rslt = tmp_tree[pos]
        if len(set_rslt) == 0:
            return False, -1
        for rcd in set_rslt:
            start_pos = rcd[0]
            return True, start_pos
        return False, -1

    def is_within_exon_region_interval_tree(self, chrm, pos):
        if chrm not in self.interval_tree_exon:
            return False
        tmp_tree = self.interval_tree_exon[chrm]
        set_rslt = tmp_tree[pos]
        if len(set_rslt) == 0:
            return False
        else:
            return True

    # 保存count矩阵
    def new_matrix(self,file_name):
        fa = open(self.out_path +file_name+".count.txt", 'w')
        fa.write('\t')
        for sub in self.TE_counter:
            fa.write(str(sub) + '\t')
        fa.write('sum_exon\tsum_map\n')
        fa.close()

    def clean_te_counter(self):
        for sub in self.TE_counter:
            self.TE_counter[sub] = 0

    # 判断read是否为重复序列，且所属family、subtype
    def te_judge(self, chrm, start):
        do = self.is_within_repeat_region_interval_tree(chrm, start)
        if do[0]:
            fam = self.te_annotation[chrm][do[1]][0][2]
            sub = self.te_annotation[chrm][do[1]][0][1]
            return 1, [fam, sub]
        return 0, [None, None]

    # 检验是否在exon上面
    def exon_judge(self, mapping_quality, exonmq, chrm, start):
        if mapping_quality >= exonmq:  # mapping quality在15之上
            on_exon = self.is_within_exon_region_interval_tree(chrm, start)
            if on_exon:
                return True
            else:
                return False

    def get_barcode_index(self, fileds, cell_barcode):
        for i in range(len(fileds)):
            if fileds[i][:2] == cell_barcode:
                return i

    def get_repeat_num(self, file, mp, temq, exonmq):   # process single sequencing file
        self.clean_te_counter() 

        if file[-1] == '/':
            file = file.rstrip('/')
        file_name = (file.split('/'))[-1]
        # file_format = file_name.split('.')[-1]  # 得到后缀，判断是什么格式的文件
        print('[%s] Start processing the BAM/SAM file...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        sum_map = 0
        sum_exon = 0
        
        bam_file = pysam.AlignmentFile(file, "rb")
        for i in bam_file:
            if i.is_duplicate or i.is_unmapped or i.is_supplementary:
                continue
            sum_map += 1
            if self.getM(i.cigarstring) >= i.query_length * mp and i.mapping_quality >= temq:
                chrm = i.reference_name
                if 'chr' not in chrm:
                    chrm = 'chr' + chrm
                start = i.pos + 1
                te_result = self.te_judge(chrm, start)
                if te_result[0]:
                    self.TE_counter[te_result[1][1]] += 1
                if self.exon_judge(i.mapping_quality, exonmq, chrm, start):
                    sum_exon += 1
        bam_file.close()

        # ---输出文件的位置---
        fa = open(self.out_path + "all.count.txt", 'a')
        fa.write(file_name + '\t')  # first at each row
        for sub in self.TE_counter:
            fa.write(str(self.TE_counter[sub]) + '\t')
        fa.write(str(sum_exon) + '\t' + str(sum_map) + '\n')
        fa.close()
        print('[%s] Results were written.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def get_repeat_num_10x_bc_sorted(self, file, mp, temq, exonmq, cell_barcode):
        if file[-1] == '/':
            file = file.rstrip('/')
        file_name = (file.split('/'))[-1]
        file_format = file_name.split('.')[-1]
        print('[%s] Start processing the BAM/SAM file...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

        self.new_matrix(file_name)

        fa = open(self.out_path + file_name+ ".count.txt", 'a')
        barcode_dic = {}
        sum_map = 0
        sum_exon = 0

        bam_file = pysam.AlignmentFile(file, "rb")  # 读取文件
        umi_set = set()
        last_barcode = ''
        for i in bam_file:  # 迭代
            try:
                barcode = i.get_tag('CB').split('-')[0]  # 获得barcode
                umi = i.get_tag('UB')
            except KeyError:
                continue

            if self.set_f_bc!=set():
                if barcode not in self.set_f_bc:
                    continue

            if barcode not in barcode_dic:  # 检查是否是新的barcode，若是新的，要写入数据还要重置
                barcode_dic[barcode] = 0

                if last_barcode != '' and len(umi_set) > 500:  # 如果不是第一次读取，则输入上一次TE_counter的计数结果
                    fa.write(last_barcode + '\t')
                    for sub in self.TE_counter:
                        fa.write(str(self.TE_counter[sub]) + '\t')
                    fa.write(str(sum_exon) + '\t' + str(sum_map) + '\n')

                if last_barcode != '':
                    for sub in self.TE_counter:
                        self.TE_counter[sub] = 0
                    sum_map = 0
                    sum_exon = 0
                    umi_set = set()

                last_barcode = barcode

            umi_set.add(umi)

            if i.is_duplicate or i.is_unmapped or i.is_supplementary:  # 判断是否是重复、未比对上、次要匹配
                continue
            sum_map += 1
            if self.getM(i.cigarstring) >= i.query_length * mp and i.mapping_quality >= temq:  # 比较质量和cigar
                chrm = i.reference_name  # 通过质量比较则获取chrm和start position
                if 'chr' not in chrm:
                    chrm='chr'+chrm
                start = i.pos + 1
                te_result = self.te_judge(chrm, start)  # te_judge是去间隔树看，返回是否是重复序列，和他的family和类别
                if te_result[0]:
                    self.TE_counter[te_result[1][1]] += 1
                if self.exon_judge(i.mapping_quality, exonmq, chrm, start):  # 比一下read是否在exon上
                    sum_exon += 1
        bam_file.close()

        fa.write(last_barcode + '\t')
        for sub in self.TE_counter:  # 最后一个cell的barcode没法在上面写入，要额外写进去
            fa.write(str(self.TE_counter[sub]) + '\t')
        fa.write(str(sum_exon) + '\t' + str(sum_map) + '\n')
        fa.close()
        print('[%s] Results were written.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def getM(self, str):
        regex = re.compile('\d{1,3}M')
        data = regex.findall(str)
        sum_M = 0
        for i in range(len(data)):
            sum_M += int(data[i].replace('M', ''))
        return sum_M

    def get_all_barcode(self, file):
        all_barcode = set()
        bam_file = pysam.AlignmentFile(file, "rb")  # 读取文件
        for i in bam_file:
            barcode = i.get_tag('CB').split('-')[0]  # 获得barcode
            all_barcode.add(barcode)
        return all_barcode

    def get_repeat_num_10x_df(self, file, mp, temq, exonmq, cell_barcode):
        cell_barcode_index = -1

        if file[-1] == '/':
            file = file.rstrip('/')
        file_name = (file.split('/'))[-1]
        file_format = file_name.split('.')[-1]
        print('Start processing the', file_format, 'file...')

        if file_format == 'sam':
            # ---输入sam文件---
            f = open(file, 'r')
            while True:
                line = f.readline().rstrip()
                if len(line) == 0:  # 结束while循环
                    break
                elif line[0] == '@':  # 以SRR开头的是序列匹配情况
                    continue
                fileds = line.split('\t')

                if cell_barcode_index == -1:
                    cell_barcode_index = self.get_barcode_index(fileds, cell_barcode)

                flag = int(fileds[1])
                barcode = fileds[cell_barcode_index].split(':')[2].split('-')[0]  # 拆分barcode CB:Z:AAAGATGAGGTCGGAT-1
                te_matrix.setdefault(barcode, [copy.deepcopy(self.TE_counter), {'sum_map': 0, 'sum_exon': 0}])
                if (flag & 4 != 0) or (flag & 1024 != 0) or (
                        flag & 256 != 0):  # 4为unmapped;1024为duplicate;256为次要匹配
                    continue
                te_matrix[barcode][1]['sum_map'] += 1
                length = len(fileds[9])
                if (self.getM(fileds[5]) >= length * mp) and (int(fileds[4]) >= temq):  # l[4]是mapping quality
                    chrm = fileds[2]
                    start = int(fileds[3])
                    te_result = self.te_judge(chrm, start)
                    if te_result[0]:
                        te_matrix[barcode][te_result[1][0]][te_result[1][1]] += 1
                    if self.exon_judge(int(fileds[4]), exonmq, chrm, start):
                        te_matrix[barcode][1]['sum_exon'] += 1
            f.close()
        elif file_format == 'bam' or file_format == 'cram':
            barcode_dic = {}
            bam_file = pysam.AlignmentFile(file, "rb")  # 读取文件
            for i in bam_file:  # 迭代
                barcode = i.get_tag(cell_barcode).split('-')[0]  # 获得barcode

                if i.is_duplicate or i.is_unmapped or i.is_supplementary:  # 判断是否是重复、未比对上、次要匹配
                    continue

                if barcode not in barcode_dic:
                    self.df.loc[barcode] = 0
                    barcode_dic[barcode] = 0

                self.df.loc[barcode, 'sum_map'] += 1

                if self.getM(i.cigarstring) >= i.query_length * mp and i.mapping_quality >= temq:  # 比较质量和cigar
                    chrm = i.reference_name  # 通过质量比较则获取chrm和start position
                    start = i.pos + 1
                    te_result = self.te_judge(chrm, start)  # te_judge是去间隔树看，返回是否是重复序列，和他的family和类别
                    if te_result[0]:
                        self.df.loc[barcode, te_result[1][1]] += 1
                    if self.exon_judge(i.mapping_quality, exonmq, chrm, start):  # 比一下read是否在exon上
                        self.df.loc[barcode, 'sum_exon'] += 1
            bam_file.close()
            print('Write the result')

        # ---输出文件的位置---
        #         fa = open(self.out_path + "count.txt", 'a')
        self.df.to_csv(self.out_path + 'count.csv')

    #         fa.close()
