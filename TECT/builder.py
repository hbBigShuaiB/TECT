import pickle

from TECT.intervaltree import IntervalTree
from TECT.tools import *
from datetime import datetime


class Builder:
    def __init__(self, rmsk, exon_gff, out):  # 使用rmsk，gff进行初始化；包含out路径

        self.out_path = out  # 输出的路径

        self.f_te_annotation = rmsk  # repeat masker注释文件的位置
        self.f_exon_annotation = exon_gff

        self.te_annotation = {}  # 每条染色体->{每个基因区域的起始位置->(end_pos,b_rc,sub_type,family,csn_start,csn_end)}
        self.interval_tree_te = {}  # 构建间隔树，字典key为chr，value里存树

        self.exon_annotation = {}  # 每条染色体->{}
        self.interval_tree_exon = {}

        self.TE_counter = {}  # 记录TE数量



        self.load_rmsk_annotation()  # load repeat注释文件
        self.index_te_interval_tree()  # 将repeat所在的位置写成间隔树保存
        self.load_exon_annotation()  # load exon注释文件
        self.index_exon_interval_tree()  # 将exon所在位置写成间隔树
        self.save_to_index_file()  # 保存为idx文件

    # 7288   1.3  0.1  0.1  chr21     32213340 32214177 (15915718) +  L1HS LINE/L1 5312 6155    (0) 4482983
    # 4741   0.6  0.0  0.0  chr21     17916705 17917238 (30212657) C  L1HS LINE/L1 (0) 6155   5622 4462908
    def load_rmsk_annotation(self):
        print('[%s] Start the TECT_build module.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        print('')
        print('[%s] Loading the TE annotation file...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        with open(self.f_te_annotation) as f_te:  # 这个是从repeatmasker上下载的重复序列所在染色体位置的文件
            for line in f_te:
                if len(line.rstrip()) < 1:
                    continue
                fields = line.split()
                if len(fields) < 15:
                    continue
                if fields[0].rstrip() == "SW" or fields[0].rstrip() == "score":
                    continue
                sub_type = fields[9]  # 子类型
                family = fields[10]  # 家族
                if family == 'Simple_repeat' or family.split('/')[0] == 'Satellite' or family == 'Low_complexity':
                    continue
                chrm = fields[4]
                start_pos = int(fields[5])
                end_pos = int(fields[6])

                # b_rc = False
                # if fields[8] == "C":  # 反链
                #     b_rc = True

                # 210905改为仅记录subtype
                # if family not in self.TE_counter:
                #     self.TE_counter[family] = {}  # 若无family，TEcounter设定family为{}
                if sub_type not in self.TE_counter:
                    self.TE_counter[sub_type] = 0

                # if b_rc:
                #     csn_start = int(fields[13])
                #     csn_end = int(fields[12])
                # else:
                #     csn_start = int(fields[11])
                #     csn_end = int(fields[12])

                if chrm not in self.te_annotation:
                    self.te_annotation[chrm] = {}
                if start_pos in self.te_annotation[chrm]:  # this is not allowed!!!!
                    print("Position {0}:{1} has more than  1 annotation!!!".format(chrm, start_pos))

                if start_pos not in self.te_annotation[chrm]:
                    self.te_annotation[chrm][start_pos] = []
                self.te_annotation[chrm][start_pos].append((end_pos, sub_type, family))
        print('[%s] TE annotation file had been loaded.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def index_te_interval_tree(self):
        print('[%s] Building the TE interval trees...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        for chrm in self.te_annotation:
            interval_tree = IntervalTree()
            for pos in self.te_annotation[chrm]:
                end_pos = self.te_annotation[chrm][pos][0][0]
                interval_tree.addi(pos, end_pos)
            self.interval_tree_te[chrm] = interval_tree
        print('[%s] The TE interval trees had been built.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def load_exon_annotation(self):
        print('[%s] Loading the exon annotation file...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        with open(self.f_exon_annotation) as f_exon:
            for line in f_exon:
                if line[0] == '#':
                    continue
                fields = line.split()
                chrm = fields[0]
                if 'chr' not in chrm:
                    chrm='chr' + chrm
                region_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                if region_type == 'exon' and start != end:
                    if chrm not in self.exon_annotation:
                        self.exon_annotation[chrm] = []
                    self.exon_annotation[chrm].append((start, end))
        print('[%s] The exon annotation file had been loaded.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def index_exon_interval_tree(self):
        print('[%s] Start Building the exon interval trees...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        for chrm in self.exon_annotation:
            interval_tree = IntervalTree()
            for pair in self.exon_annotation[chrm]:
                start = pair[0]
                end = pair[1]
                interval_tree.addi(start, end)
            self.interval_tree_exon[chrm] = interval_tree
        print('[%s] The exon interval trees had been built.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    def save_to_index_file(self):
        print('[%s] Start generating the index file...'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        if self.out_path[-1] != '/':
            self.out_path += '/'

        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)

        self.out_path = self.out_path + 'index.idx'
        idx_file = {'te': self.te_annotation, 'te_tree': self.interval_tree_te,
                    'exon_tree': self.interval_tree_exon, 'counter': self.TE_counter}
        with open(self.out_path, 'wb') as f_index:
            pickle.dump(idx_file, f_index)
        print('[%s] The index file was built successfully.'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
