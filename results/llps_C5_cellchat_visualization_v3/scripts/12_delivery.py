from pathlib import Path
import csv,hashlib,zipfile,json
O=Path(__file__).resolve().parents[1];D=O/'delivery';D.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
(D/'UPLOAD_STATUS.md').write_text('''# 上传状态

2026-11-26 已通过 GitHub connector 尝试向 `Kimsx05/transfer_for_LLPS` 的 `results/llps_C5_cellchat_visualization_v3/README.md` 新建文件。

结果：失败，HTTP 403，`Resource not accessible by integration`。没有成功上传文件或提交。没有尝试读取本机缓存凭据或绕过仓库权限。

本地完整轻量交付包：`llps_C5_cellchat_visualization_v3_lightweight.zip`。包含图片、绘图数据、统计表、报告、脚本及历史图归档；没有 Seurat、CellChat 对象或完整表达矩阵。解压后打开 FIGURE_GALLERY.html 可离线浏览。

压缩包内目录前缀与目标仓库一致：`results/llps_C5_cellchat_visualization_v3/`。output_manifest.csv 列出内容及 SHA256 (不包含清单自身及压缩包，以避免递归)。压缩包 SHA256 另见 package_sha256.txt。
''')
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name not in ['output_manifest.csv','llps_C5_cellchat_visualization_v3_lightweight.zip','package_sha256.txt'])
assert not any(p.suffix.lower() in ['.rds','.qs','.qs2','.h5ad','.h5seurat'] for p in files)
with (O/'output_manifest.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['relative_path','absolute_path','bytes','sha256'])
 for p in files:w.writerow([p.relative_to(O),p,p.stat().st_size,sha(p)])
zpath=D/'llps_C5_cellchat_visualization_v3_lightweight.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files+[O/'output_manifest.csv']:z.write(p,'results/llps_C5_cellchat_visualization_v3/'+str(p.relative_to(O)))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
(D/'package_sha256.txt').write_text(sha(zpath)+'  '+zpath.name+'\n')
print('PACKAGED',len(files)+1,'files',zpath.stat().st_size,'bytes',zpath)
