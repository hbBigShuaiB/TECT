# 对count matrix，选择一种聚类方法，按不同的细胞类型进行聚类

import subprocess
import pandas as pd
import scanpy as sc
import numpy as np
import os
from matplotlib import pyplot as plt

from TECT.tools import *

class Cluster:
    def __init__(self, TEMatrix,geneMatrix,out_path,max_count,max_gene,mt_percent,top_gene,gene_resolution,te_resolution):
        self.te_matrix=TEMatrix
        self.gene_matrix=geneMatrix

        self.gene=None
        self.te=None

        self.max_count=max_count
        self.max_gene=max_gene
        self.mt_percent=mt_percent
        self.top_gene=top_gene
        self.gene_resolution=gene_resolution
        self.te_resolution=te_resolution
        
        self.out_path=out_path
        if self.out_path[-1] != '/':
            self.out_path += '/'
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)

        # self.command = 'Rscript'
        # self.path2script = ''
        # self.cmd = ''

        self.setting()
        self.gene_process()
        self.te_process()

    def setting(self):
        sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=150,dpi_save=300,facecolor='white')
        sc.settings.figdir=self.out_path

    def read_file(self,f_path):
        if os.path.isdir(f_path):
            adata=sc.read_10x_mtx(f_path,'gene_symbols',make_unique=True,cache=True)
        else:
            if f_path[-1] == '/':
                f_path = f_path.rstrip('/')
            file_name = (f_path.split('/'))[-1]
            file_format = file_name.split('.')[-1]  # 得到后缀，判断是什么格式的文件
            if file_format=='txt':
                adata=sc.read_text(f_path)
        return adata

    def do_pca(self,adata):
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        return adata

    def gene_process(self):
        gene=self.read_file(self.gene_matrix)

        if '-' in gene.obs.index[0]:
            # remove sample lable
            newIndex=gene.obs.index.tolist()
            for i in range(len(newIndex)):
                newIndex[i]=newIndex[i].split('-')[0]
            gene.obs.index=newIndex

        # sc.pp.filter_cells(gene, min_genes=200)
        sc.pp.filter_genes(gene, min_cells=3)

        # Mitochondrial content
        gene.var['mt'] = gene.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(gene, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        ##Another to save figures
        # if not os.path.exists(self.out_path+'preprocess'):
        #     os.makedirs(self.out_path+'preprocess')

        # with plt.rc_context():
        #     # violin
        #     sc.pl.violin(gene, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        #         jitter=0.4, multi_panel=True,show=False)
        #     plt.savefig(self.out_path+"preprocess/01.violin.png",dpi=150, bbox_inches="tight")

        ##1.violin & 2.scatter
        # show the gene expression, used to adjust the cell filtering threshold
        sc.pl.violin(gene, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=0.4, multi_panel=True,save='.png')
        sc.pl.scatter(gene, x='total_counts', y='pct_counts_mt',save='.png')
        sc.pl.scatter(gene, x='total_counts', y='n_genes_by_counts',save='2.png')
        # filter
        sc.pp.filter_cells(gene,max_counts=self.max_count)
        sc.pp.filter_cells(gene,max_genes=self.max_gene)
        # gene = gene[gene.obs.n_genes_by_counts < 2500, :]
        gene = gene[gene.obs.pct_counts_mt < self.mt_percent, :]

        # library-size correct,so that counts become comparable among cells
        sc.pp.normalize_total(gene,target_sum=1e6)
        sc.pp.log1p(gene)

        # highly variable genes
        sc.pp.highly_variable_genes(gene, n_top_genes=self.top_gene)
        gene = gene[:, gene.var.highly_variable]

        # sc.pp.regress_out(gene, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(gene, max_value=10)

        # PCA
        gene=self.do_pca(gene)

        sc.tl.leiden(gene,resolution=self.gene_resolution)
        print(gene.obs)
        sc.pl.umap(gene, color=['leiden'],title='resolution='+str(self.gene_resolution),show=False,
            save='_gene.png')
        
        sc.tl.rank_genes_groups(gene, 'leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(gene, n_genes=25, sharey=False,show=False,save='_before.png')

        # # leiden
        # with plt.rc_context():
        #     sc.tl.leiden(gene,resolution=0.1)
        #     sc.pl.umap(gene, color=['leiden'],title='resolution=0.1',show=False)
        #     plt.savefig(self.out_path+'clust/01.UMAP.png', bbox_inches="tight")

        #     sc.tl.rank_genes_groups(gene, 'leiden', method='wilcoxon')
        #     sc.pl.rank_genes_groups(gene, n_genes=25, sharey=False,show=False)
        #     plt.savefig(self.out_path+'clust/02.markerTE.png',bbox_inches='tight')
        self.gene=gene

    def te_process(self):
        gene=self.gene

        te=self.read_file(self.te_matrix)

        # drop and store sum count
        te.obs['sum_map']=te.X[:,-1].tolist()
        te.obs['sum_exon']=te.X[:,-2].tolist()
        te=te[:,0:te.n_vars-2]

        te.obs['n_counts'] = te.X.sum(1)
        te.obs['n_genes'] = (te.X > 0).sum(1)
        
        te=te[gene.obs.index]
        sc.pp.filter_genes(te, min_cells=3)

        cluNo=gene.obs.leiden.cat.categories.to_list()
        newNo=len(cluNo)
        gene.obs['te_leiden']=gene.obs['leiden']

        # clusters by gene result
        for i in cluNo:
            cluster=te[gene.obs.leiden==i]
            sc.pp.normalize_total(cluster) 
            sc.pp.log1p(cluster)
            sc.pp.highly_variable_genes(cluster,min_disp=0,max_mean=float('inf'),min_mean=0)

            cluster = cluster[:, cluster.var.highly_variable]
            sc.pp.scale(cluster, max_value=10)

            cluster=self.do_pca(cluster)

            sc.tl.leiden(cluster,resolution=self.te_resolution)
            each_clust_leiden= cluster.obs.leiden

            if len(each_clust_leiden.cat.categories)>1:
                gene.obs.te_leiden=gene.obs.te_leiden.cat.remove_categories(i)
                gene.obs.te_leiden=gene.obs.te_leiden.cat.add_categories([str(newNo+int(g)) for g in each_clust_leiden.cat.categories])
                each_clust_leiden.cat.categories=[str(newNo+int(g)) for g in each_clust_leiden.cat.categories]
                newNo+=len(each_clust_leiden.cat.categories)
                # gene.obs['te_leiden']=cluster.obs['leiden']
                for te_clust_index in cluster.obs.index:
                    gene.obs.at[te_clust_index,'te_leiden']=cluster.obs.at[te_clust_index,'leiden']

                sc.tl.rank_genes_groups(cluster, 'leiden',method='wilcoxon')
                sc.pl.rank_genes_groups(cluster, n_genes=25, sharey=False,save='_te_cluster'+i+'.png')
            
                sc.pl.umap(cluster, color=['leiden'],title='cluster the cluster' +i + ' with te expresion',save='_te_cluster'+i+'.png')
            
        sc.pl.umap(gene, color=['te_leiden'],title='re-project with the te result',save='_all_reproject.png')
        sc.tl.rank_genes_groups(gene, 'te_leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(gene, n_genes=25, sharey=False,show=False,save='_after.png')

        self.gene=gene
        self.te=te

    def km(self):
        self.path2script = 'TECT/kmeans.R'

        self.cmd = [self.command, self.path2script] + self.args
        subprocess.run(self.cmd, universal_newlines=True)

    def km_consensus(self):
        mkdir('result')
        mkdir('result/km')
        self.path2script = 'TECT/km_consensus.R'

        self.cmd = [self.command, self.path2script] + self.args
        subprocess.run(self.cmd, universal_newlines=True)

    def hc_consensus(self):
        mkdir('result')
        mkdir('result/hc')
        self.path2script = 'TECT/hc_consensus.R'

        self.cmd = [self.command, self.path2script] + self.args
        subprocess.run(self.cmd, universal_newlines=True)
