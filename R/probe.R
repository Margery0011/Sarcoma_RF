probe = read.table(file = 'GPL13534-11288.txt',
                   sep = '\t',
                   quote = '',
                   comment.char = '#', # 过滤掉'GPL13534-11288.txt'文件中以‘#’开头的注释信息
                   header = T,
                   fill = T,           #  如果文件中某行的数据少于其他行，则自动添加空白域。     
                   stringsAsFactors = F) # 字符串不改为因子

ids = probe[probe$Symbol != '',
            c(1,13)] # 提取探针和geneID

a=read.table(file = "GSE140686-GPL13534_series_matrix.txt.gz")
