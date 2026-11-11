from env import *
import html,zipfile
E=read('07_summary/evidence_matrix.tsv.gz');mf=read('01_signatures/signature_manifest.tsv');v=json.loads((B/'00_plan_audit/validation.json').read_text());fig=read('07_summary/figure_index.tsv');c=E[(E.module=='clinical')&(E.signature=='MP10')&E.main.fillna(False)&E.adjustment.eq('unadjusted')&E.P.notna()];cm=E[(E.module=='cell_models')&(E.signature=='MP10')&E.effect_type.eq('point_biserial_r')];ep=E[(E.module=='epi')&(E.dataset=='InhouseData')&E.subset.eq('Tumor')];drug=E[(E.module=='measured_drug')&(E.signature=='MP10')&E.method.eq('spearman')&E.endpoint.eq('LN_IC50')];night=json.loads((B/'logs/night_status.json').read_text()) if (B/'logs/night_status.json').exists() else {'state':'not_launched'}
parts=[]
def p(x):parts.append('<p>'+x+'</p>')
def h(x):parts.append('<h2>'+x+'</h2>')
def table(d,columns=None):
 if columns:d=d[[x for x in columns if x in d]]
 parts.append('<div class="table">'+d.to_html(index=False,escape=True,float_format=lambda x:f'{x:.4g}',na_rep='不可估')+'</div>')
def link(label,path):return '<a href="../'+path+'">'+html.escape(label)+'</a>'
def imagefig(stem):parts.append(f'<figure><a href="figures/{stem}.pdf"><img src="figures/{stem}.png" loading="lazy"></a><figcaption>{html.escape(stem)} (PDF/300 dpi PNG；完整解释见图索引)</figcaption></figure>')
parts.append('<h1>MP10 耐药关联全量扩展验证 v2</h1>');p('生成时间：'+datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()+'。这是看过 v1 后开展并冻结的探索性扩展，不能称为事前预注册；不以显著阳性为验收条件。')
p('<strong>当前证据不支持把 MP10 定义为跨背景通用的顺铂耐药指标。</strong> 现有结果更符合背景和状态依赖：部分获得性模型上升，另一些下降；Epi 的两个 Tsunoda 方向并不同时支持耐药假设；实测顺铂药敏和临床主要终点未出现稳定支持。所有反向和不显著结果保留。')
if not v['official_full_scoring_complete']:
 p('<strong>尚未全部结束：</strong>average-tie Epi、真实周期调整、pseudobulk、细胞模型、实测 GDSC、可评估临床及真实配对分支已执行。官方 AUCell/UCell 全量方法敏感性仍为 DEFERRED_RESOURCE；256 细胞真实试运行成功，排名约11.5秒，估算单进程约93分钟，按项目规定留在23:00–07:00 Asia/Shanghai。当前本地调度状态：'+html.escape(str(night))+ '。不能将排队写成完成。夜间脚本结束后自动重算方法比较、FDR、图及本报告和审阅包。')
else:p('官方 AUCell/UCell 全量方法敏感性亦已完成；按用户追加授权改为立即运行，统计定义未变。逐细胞ID、完整性、有限分数及周期重算一致性均已核验；具体方法间差异见 Epi 全量表和图。')
h('执行范围与独立单位')
p(f"登记 {v['signatures_registered']} 个评分定义/敏感性版本，其中 {v['signatures_members_ready']} 个获得成员或按冻结规则构建；其余源定义阻塞。版本包含相同成员的敏感性别名，不应视为独立证据。Epi 使用123,699细胞、57样本、5数据集；只核实 Inhouse 12样本对应11患者，其中达到肿瘤Epi最低50细胞条件的患者为 {int(ep[ep.method.eq('average_tie')].n.max())} 位。其余数据集照常评分、描述，但不虚构患者推断。")
p(f"细胞模型覆盖 {v['studies']} 项研究、{v['paired_backgrounds']} 种有亲本比较的细胞背景、{v['contrasts']} 个状态/阶段比较，按阶段折叠后为 {v['study_background_drug_process_units']} 个研究内背景—耐药过程单元；这不等于独立诱导谱系。原始RNA测定 {v['cell_model_assays']} 次，其中兼容的基线/DMSO测定 {v['baseline_model_assays']} 次；RT112只有亲本表达。GRC1/2的P0同名且超过99.98%共有基因数值一致，按依赖控制处理；物理RNA是否同份仍需原始对应，不能增加独立背景数。")
p(f"有至少一个主要临床模型可评估的患者ID，按队列内去重后合计 {v['clinical_unique_IDs_with_MP10_evaluable_primary']} 个，分布于6项研究、7个治疗人群；Taber NAC/晚期分开，同一患者不跨这两个主分支重复。该总数不能证明不同研究间完全没有隐匿患者重叠。其余队列完成评分/QC，疗效或配对标签缺失则单独阻断。")
table(read('05_clinical/cohort_counts.tsv'),['dataset','n_assays','n_patient_IDs','branch','endpoint','n_eligible','n','events','responders','status'])
h('冻结方案、评分和 v1 技术修复')
p(link('冻结配置','00_plan_audit/frozen_config.json')+'；'+link('全部计划检验','00_plan_audit/hypothesis_table.tsv')+'；'+link('覆盖表','00_plan_audit/coverage.tsv')+'；'+link('原生成员/权重','01_signatures/members.tsv')+'；'+link('每队列实际评分版本','01_signatures/implementation_by_dataset.tsv')+'。MP10 Top20固定，主分析连续评分，无NMF重跑、最佳切点、表达预测IC50或生存替代。80%覆盖是本轮操作阈值，不是公认生物学界值。未测基因与零表达严格区分。')
p('counts采用TMM及logCPM(prior.count=1)，线性TPM/CPM/RPKM用log2(x+1)，作者已转换矩阵不再log。完整背景排名；逐集合ssGSEA先算未归一化值，再除以同一评分人群的本集合极差，逐队列核对等于v1单集合normalize=TRUE，避免多集合共同归一化改变MP10单位。native逐基因z只使用冻结参考人群；score_z也不因终点或敏感性子集重算。所有partial/adapted有单独score_name；六基因比值必须全覆盖且有明确非负线性尺度。SCAN、未知线性基准和不清楚offset的作者logCPM不计算比值。')
table(read('00_plan_audit/bugfix_audit.tsv'));p('v1保护核验全部通过：'+link('逐文件SHA256对照','00_plan_audit/v1_protection_check.tsv')+'。本地不是Git仓库，以实际文件hash固定版本；参考commit已记录，未声称本地等于远程commit。')
table(read('00_plan_audit/v1_v2_model_before_after.tsv'))
p('GSE212810 的 v2 OR 改变主要来自事前冻结的小样本Firth规则；同一新旧分数使用同一MLE后的对照可将评分修复与估计方法区分。未按P值选择估计器。Epi原MP10 average-tie分数复现最大绝对误差低于1e-12；它不冒称官方AUCell。')
p('冻结后另有一项明确记录的实现修订：GSE285648原生基因z参考人群收窄为15个DMSO基线样本，PARP处理的21个样本继续单列评分/QC。原先比较已排除PARP，但基因z参考曾混入，已保存修复前后分数；该修正不依据P值，MP10单样本ssGSEA及其主要模型方向不因此改变。'+link('修订记录','00_plan_audit/amendments/GSE285648_DMSO_reference.json'))
h('Epi：转录关联与调整后的衰减')
table(ep[(ep.signature.isin(['TSUNODA_UP','TSUNODA_DN','TSUNODA_UP_no_INSIG1','TSUNODA_UP_no_C1S','TSUNODA_DN_no_C1S','TSUNODA_UP_no_INSIG1_C1S']))&(ep.method=='average_tie')],['signature','adjustment','n','effect','CI_low','CI_high','P','FDR','direction'])
p('每样本内连续Spearman或rank-residual相关；同一患者a/b样本先等权Fisher-z合并，再患者等权。CI为整患者bootstrap的平均Fisher效应区间；P是双侧精确方向检验，原假设为患者效应为正的概率等于1/2，因此CI与P针对不同概念，可能出现均值CI不含0但方向检验不显著。绝不采用123,699个细胞作为独立重复的P值。')
p('UP为正但DN也为正，后者与预期耐药方向相反；深度/周期调整后两者均衰减。衰减可由深度、共同细胞状态及真实生物学共同造成，不能自动说原结果是假象。该模块只支持转录关联，不证明某个Epi细胞已耐药。原始Type分层、正常来源、既有CNV子集、nFeature替代、leave-one-patient/dataset均见全表。只有一个数据集可做可靠患者推断，leave-one-dataset移除它后无法估计，这是身份资料限制。')
imagefig('A_Tsunoda_individual_patients');imagefig('A_all_datasets_adjustments')
table(ep[(ep.method=='pseudobulk_native')&ep.signature.isin(['CisSig19','Xylinas9','NPEPPS_VRAC_ratio','Tanaka49','CrM63_weighted_z','CLSPN_derived_program'])],['signature','n','effect','CI_low','CI_high','P','FDR'])
p('pseudobulk先按sample汇总原始counts，再按dataset独立归一化；它不等于患者平均单细胞分数。100个固定随机集合按 pooled pseudobulk 表达和检出率五分位匹配，作为描述性背景对照，不能说是单细胞检测率匹配或已证实耐药特异性。'+link('随机集合与结果','02_epi/random_control_associations.tsv'))
if v['official_full_scoring_complete']:
 h('官方 AUCell/UCell 全量方法敏感性')
 concordance=read('02_epi/method_score_concordance.tsv');concordance=concordance[(concordance.signature=='MP10')&(concordance.n_cells>=50)]
 table(concordance.groupby(['method_a','method_b']).Spearman.agg(n_samples='count',median='median',minimum='min',maximum='max').reset_index())
 p('上表是达到50细胞条件的样本内MP10分数Spearman一致性，属于描述性方法比较；不是患者独立重复的效应检验。')
 table(ep[(ep.signature.isin(['TSUNODA_UP','TSUNODA_DN']))&ep.method.isin(['official_AUCell','UCell'])&ep.adjustment.isin(['raw','depth_cycle'])],['signature','method','adjustment','n','effect','CI_low','CI_high','P','FDR','direction'])
 p(link('全量完成及既有结果保护核验','00_plan_audit/full_completion_check.json')+'；'+link('实际4 CPU上限与环境验收','00_plan_audit/full_execution_resource_verification.json'))
 p('官方AUCell保留了UP正相关及DN反向正相关的模式，UCell两项效应更弱；深度/周期调整后，三种方法的UP、DN均未通过各自家族FDR。官方方法结果没有把原先的相反方向或不显著结果转成一致耐药证据。')
 p(link('逐样本全部集合的方法一致性','02_epi/method_score_concordance.tsv')+'；'+link('官方全量验收','00_plan_audit/official_full_QC.json')+'。方法一致性为描述性分析，未把细胞数当独立生物学重复；不同方法的分数尺度不直接互换。')
h('细胞模型：一致与相反方向均保留')
table(cm,['dataset','background','state','n_parental','n_resistant','delta','delta_low','delta_high','effect','CI_low','CI_high','g','g_low','g_high','P','FDR'])
p('可复现的正向包括GSE231835 T24R2，以及GSE171537的5637 CisR；反向包括UM-UC-3 CR、RT4R、T24LTT/J82LTT等。TCCSUP、NTUB1等点估计虽向上，也必须连同区间和FDR解读。T24在不同实验室/衍生状态可相反，不能只挑一种当验证成功。RT4R取样带药与亲本PBS的差异包含急性暴露混杂；UM-UC-3维持用药、洗脱未明确；T24R2原文明确维持6.6μM顺铂，但RNA取样前是否洗脱没有独立说明；GSE171537则明确多次无药传代。Base47/驱动基因无逐株可靠资料处保持未知，不把空白当WT。')
p('吉西他滨阶段不是“天”：T24 P3/P7下降、P15回升；5637两种建立过程中多为下降或不确定，且P0近乎相同并按依赖处理。该模式更接近状态依赖，不能自动外推顺铂。point-biserial r及Hedges g量化同背景RNA重复间的组分离，不是独立耐药演化重复；极小组内方差可产生极大g。GSE58624为1+1，只给差值，不制造r/P。各研究归一化分数差值不直接跨研究平均。')
table(read('03_cell_models/GSE171537_background_summary.tsv').query('signature=="MP10"'));p('原五背景CisR两升三降，精确符号检验未形成总体支持；五个效应即便全同向，双侧符号检验最小P也只能0.0625。')
imagefig('B_all_MP10_model_effects');imagefig('B_all_signature_model_r');imagefig('B_MP10_gene_changes')
h('实测药敏与耐药倍数')
table(drug,['dataset','n','effect','CI_low','CI_high','P','FDR']);p('GDSC2主分析17个匹配稳定SIDM背景，GDSC1为18个；两者高度重叠，不合计成35独立背景。MP10与实测LN_IC50无稳定关联。CisSig训练来源包含GDSC2，其结果只算来源重叠的内部一致性。所有候选同时执行：Spearman主分支没有通过本轮FDR的关联；GDSC1补充Pearson中的Xylinas9和DNA repair有信号，但不是MP10证据，也不能择优替代主方法。AUC原始方向保留，数值更高表示同一剂量窗口下存活更强。')
table(read('04_measured_drug_response/paired_IC50_source_values.tsv'));p('原文实测表确认T24/J82的亲本与LTT IC50，以及T24R2的文中10→100μM值。每研究可匹配背景仅2或1个，不进行有意义的跨背景相关推断或跨测定无条件合并。GSE171537补表ST1已成功取得，但它明确是维持用药浓度，不是IC50，故没有误算耐药倍数。RT112有药敏表但没有耐药株RNA基线，不能伪造配对。')
imagefig('E_GDSC_MP10_measured_IC50');imagefig('E_paired_measured_IC50_points')
h('临床：真实病理与作者反应层分别报告')
table(c,['dataset','branch','endpoint','n','events','responders','effect','CI_low','CI_high','P','FDR','estimator'])
p('OR事件为nonresponse或non_pCR；每增加1 SD分数。各主要MP10模型CI跨1且家族FDR未达到0.05；GSE169455/GSE212810点估计反而偏向更好反应，GSE52219/Taber NAC良好病理点估计偏向更差反应，不能据此断言已建立临床异质性。全候选临床模型本轮也没有通过FDR的结果。少量cT、亚型、方案或周期调整分别执行；缺项/稀疏/不满足信息量规则的模型明确阻断，未堆入小样本模型。')
p('GSE169455只有125例作者NAC进入主fallback；24 induction单列。没有把76GC/46ddMVAC/3other总体构成分配到患者。GSE212810保留35顺铂病例和4GCb的独立敏感性；良好病理可评估28，严格pCR29，保留未知N和原作者标签，低分病例无选择性删除。GSE52219有6良好反应但仅4严格pCR。GSE69795严格BLCA为28：排除5UTUC、1URETHRA、3cT1、1未手术；数字pstage只用于作者分期组反应，不伪造ypTN。两个UTUC TaN0拼写的作者标签补齐有独立修订记录，未纳入严格BLCA。GSE247185作者8/5标签不当作严格pCR，Patient12免疫联合有排除敏感性，“gentamicin”原文保留。')
p('Taber公开补表恢复了96例新鲜冻存表达而不是把RNA121或总体300当样本量；21例NAC和75例晚期分开，晚期明确顺铂方案主分支67例，其他方案仍做作者全人群敏感性。原始采样逐例日期不足，因此晚期结果属于作者人群反应关联，不能声称已核实每例都在首次晚期治疗前即时取材。')
p(f"MP10结局翻转核验 {v['encoding_flip_passed']}/{v['encoding_flip_checks']} 通过，OR互为倒数、P不变。GSE212810逐患者排除诊断保留同一标准化尺度，不以它决定删样本。"+link('影响诊断','05_clinical/GSE212810_MP10_leave_one_patient.tsv'))
imagefig('D_clinical_OR_unadjusted');imagefig('D_all_clinical_signature_OR')
h('治疗后状态与尚未恢复的标签')
table(E[(E.module=='paired')&(E.signature=='MP10')],['dataset','n','effect','CI_low','CI_high','P','FDR','status','reason']);p('GSE69795按GEO原始配对字段恢复23对，MP10平均变化的区间跨0；该分支包含原研究组织来源，术后残留选择及组成变化必须保留解释。GSE124305的133个术后表达已评分，但与GSE87304的真实患者配对表未取得，不能靠数字ID相似拼出论文116对。')
imagefig('F_actual_pre_post_pairs')
blocked=read('00_plan_audit/blockers.tsv.gz');table(blocked[['module','dataset','status','reason']].drop_duplicates().query('module in ["strict_branch","signature_source"]'))
p('源定义仍缺：DTP124的原始带方向Table S6及可信人鼠一对一映射；CLSPN作者2284-UP确切成员（现有PDF是全13743基因CPM表）；Choi原分类器模型/centroids。CLSPN另有事前规则派生程序，源数据内只算内部比较。Tanaka原表的组合探针及TNFIP2/SLC38A不擅自猜修；CrM保留63个原权重，weighted mean z明确标记项目适配，归于吉西他滨演进/运动性。缺少的临床标签与需要补哪一张表详见上述状态表及完整blockers。')
p('GSE235067的41基因panel未被当作全转录组。GSE244266去_F2后得到182个基础标题组，按与结果无关的表达中位数合并，189原测定保留；这182组不是已独立核实的患者数。7组重复的全背景表达相关、方案及分期一致性已审计，但仍需作者patient–GSM对应及F2的物理重复定义。GSE87304、GSE5287、术后GSE124305没有response的项目仍完成表达与覆盖QC；未用OS/TIDE/预测IC50替代真实疗效。')
p('原始Epi记录的T/N/M逐样本审计及全部可用T_stage分层已执行；临床病理肿瘤含量缺少逐患者定量字段，因此未以表达推断纯度冒充实际病理比例。'+link('Epi分期/治疗时间审计','00_plan_audit/epi_stage_timing_audit.tsv')+'；'+link('肿瘤含量调整缺口','00_plan_audit/purity_adjustment_audit.tsv'))
p(link('GSE244266重复测定QC与身份限制','00_plan_audit/GSE244266_repeat_QC.tsv'))
h('能回答与不能回答；下一步优先级')
p('可以回答：MP10在若干明确模型/状态中发生真实转录变化，且同一细胞背景的不同衍生模型可相反；目前缺少实测药敏及可比临床终点的一致支持。不能回答：MP10是否驱动耐药、能否预测顺铂治疗获益、所有Epi是否恶性、胆固醇含量/通量是否增加。甲羟戊酸/胆固醇相关转录活性不等于脂质量，更不等于胆固醇造成耐药。')
p('建议降低“泛顺铂耐药临床标志物”的投入优先级。优先选择T24R2上升与T24LTT下降的对照模型，在一致无药取样、真实独立诱导和相同药敏测定条件下复核MP10及脂质/通量，再考虑干预与救援实验。若只在特定状态成立，则限定结论范围；TME/空间关联作为状态解释问题，与耐药因果分开。资料补充优先级为临床逐例疗效/方案/取样日期与配对表，然后是未恢复的原始signature定义；不再添加一批分数寻找阳性。')
h('完整交付、重跑与资源')
for label,path in [('完整证据矩阵','07_summary/evidence_matrix.tsv.gz'),('全部状态','00_plan_audit/analysis_status.tsv.gz'),('冻结与实际FDR数','00_plan_audit/FDR_actual_vs_plan.tsv'),('来源文件hash','00_plan_audit/source_file_manifest.tsv'),('Epi患者原始效应','02_epi/patient_effects.tsv'),('细胞模型全结果','03_cell_models/contrast_results.tsv'),('临床逐例审计','05_clinical/patient_endpoint_audit.tsv'),('每模型患者清单','05_clinical/model_patient_usage.tsv.gz'),('signature全清单','01_signatures/signature_manifest.tsv'),('重跑脚本','scripts/run_all.sh'),('全量任务入口(保留夜间默认门禁)','scripts/night_pipeline.sh')]:p(link(label,path))
p(f"实测单进程峰值RSS约 {v['max_single_process_RSS_GB']:.2f} GB；全量重算采用3个单线程R工作进程及1个调度CPU余量，整个进程树共享4 CPU亲和性硬上限。启动时发现系统Rprofile将BLAS线程改为24，已中止该启动并为本任务禁用该配置后重启；早期仅设置环境变量的线程声明不视为已验证。未整体升级R环境。目录当前约 {v['disk_bytes']/1024**3:.2f} GiB；审阅包排除原始矩阵、运行环境及大细胞分数。完整PDF/300dpi PNG位于07_summary/figures；有图例者另附panel_only与legend_only PDF。原始cache和v1均只读。")
table(fig,['figure','source','note']);h('主要原始来源')
for title,url in [('MSigDB Tsunoda','https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/TSUNODA_CISPLATIN_RESISTANCE_UP.html'),('CisSig','https://www.nature.com/articles/s41698-023-00375-y'),('Xylinas','https://pmc.ncbi.nlm.nih.gov/articles/PMC5039423/'),('NPEPPS','https://pmc.ncbi.nlm.nih.gov/articles/PMC11641020/'),('Tanaka','https://pmc.ncbi.nlm.nih.gov/articles/PMC5085609/'),('CrM','https://doi.org/10.1016/j.canlet.2024.217339'),('DTP','https://pmc.ncbi.nlm.nih.gov/articles/PMC8437243/'),('Taber公开表达和临床补表','https://www.nature.com/articles/s41467-020-18640-0'),('GDSC release8.5','https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/'),('LTT实测IC50','https://pmc.ncbi.nlm.nih.gov/articles/PMC13042512/'),('T24R2原文','https://pmc.ncbi.nlm.nih.gov/articles/PMC10653028/')]:p('<a href="'+url+'">'+title+'</a>')
for g in sorted(set(read('03_cell_models/model_identity.tsv').dataset)|set(read('05_clinical/patient_endpoint_audit.tsv').dataset.str.extract('(GSE[0-9]+)',expand=False).dropna())):p(f'<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}">{g} GEO原始记录</a>')
css='body{font:16px/1.65 Arial,"Noto Sans CJK SC",sans-serif;color:#202020;max-width:1320px;margin:35px auto;padding:0 25px;background:white}h1{font-size:30px}h2{font-size:23px;border-top:1px solid #ddd;padding-top:20px}.table{overflow:auto;margin:18px 0}table{border-collapse:collapse;font-size:13px;white-space:nowrap}th,td{border:1px solid #ddd;padding:6px 10px;text-align:left}th{background:#f3f5f7}img{max-width:100%;height:auto}a{color:#2069a1}figcaption{font-size:13px;color:#555}'
(B/'07_summary/中文结果报告.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>MP10 v2 expanded</title><style>'+css+'</style><body>'+''.join(parts)+'</body></html>')
readme=f'''MP10 resistance validation v2 expanded
中文报告：{B}/07_summary/中文结果报告.html
全量结果：{B}/07_summary/evidence_matrix.tsv.gz
审阅包：{B}/mp10_resistance_validation_v2_expanded_review.zip
重跑：bash '{B}/scripts/run_all.sh'
官方AUCell/UCell默认夜间：bash '{B}/scripts/night_pipeline.sh' (23:00–07:00 Asia/Shanghai gate)
本轮即时运行授权及资源修复见00_plan_audit/amendments/immediate_full_scoring.json；按该授权重跑可设MP10_RUN_NOW=1 MP10_EPI_WORKERS=3。
状态：官方全量评分{'完成' if v['official_full_scoring_complete'] else '尚在夜间队列；pilot已完成，不代表全量完成'}。其他各项状态见00_plan_audit/analysis_status.tsv.gz。
先阅读冻结配置及source hashes。默认复用已核验公共输入/评分检查点；RECOMPUTE=1重算bulk分数。原始Epi cache及v1必须在登记的绝对路径；小型审阅包用于独立阅读，不含大型输入矩阵。
所有阻塞均有具体文件/字段，阴性/相反方向未删除。没有OS/表达预测IC50替代真实疗效。
''';(B/'README.txt').write_text(readme)
# Preserve directory hierarchy for standalone HTML and audit reading; avoid huge raw matrices/runtime.
review=B/'mp10_resistance_validation_v2_expanded_review.zip';selected=[B/'README.txt']
for folder in ['scripts','01_signatures','03_cell_models','04_measured_drug_response','05_clinical','06_paired_states','07_summary']:
 for pth in (B/folder).rglob('*'):
  if not pth.is_file() or '__pycache__' in str(pth) or pth.suffix=='.rds':continue
  if pth.name.endswith('.png') or 'counts_input' in pth.name:continue
  if pth.stat().st_size<=8*1024**2:selected.append(pth)
for pth in (B/'00_plan_audit').glob('*'):
 if pth.is_file() and pth.stat().st_size<5*1024**2 and not ('_gene_mapping' in pth.name or '_background_qc' in pth.name):selected.append(pth)
for pth in (B/'00_plan_audit/amendments').glob('*'):
 if pth.is_file() and pth.stat().st_size<5*1024**2:selected.append(pth)
for pth in (B/'02_epi').glob('*'):
 if pth.is_file() and pth.stat().st_size<5*1024**2 and not any(x in pth.name for x in ['counts','cell_','random_matched_scores']):selected.append(pth)
for pth in (B/'logs').glob('*.log'):
 if pth.stat().st_size<500000:selected.append(pth)
# HTML package points at PDFs for scientific figures when PNGs are intentionally excluded.
with zipfile.ZipFile(review,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for pth in sorted(set(selected)):
  if pth.name=='中文结果报告.html':
   text=pth.read_text();text=re.sub(r'<img src="figures/([^\"]+)\.png" loading="lazy">',r'打开正式矢量图: \1.pdf',text);z.writestr(str(pth.relative_to(B)),text)
  else:z.write(pth,pth.relative_to(B))
(B/'07_summary/review_package.json').write_text(json.dumps(dict(path=str(review),bytes=review.stat().st_size,n_files=len(set(selected)),created=datetime.datetime.now().isoformat()),indent=2));print('REPORT',B/'07_summary/中文结果报告.html','REVIEW_BYTES',review.stat().st_size,flush=True)
