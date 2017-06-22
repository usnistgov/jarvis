from pymatgen.core.periodic_table import get_el_sp, Element
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
from pymatgen.util.plotting_utils import get_publication_quality_plot
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib_venn import venn3, venn3_circles
from matplotlib.gridspec import GridSpec
the_grid = GridSpec(3, 2) 
new=['mp-1821', 'mp-22871', 'mp-22870', 'mp-22872', 'mp-22874', 'mp-22876', 'mp-2809', 'mp-22878', 'mp-632475', 'mp-23052', 'mp-23057', 'mp-19565', 'mp-29985', 'mp-19568', 'mp-19744', 'mp-11595', 'mp-23332', 'mp-27346', 'mp-27340', 'mp-19585', 'mp-27343', 'mp-22272', 'mp-27182', 'mp-28940', 'mp-30034', 'mp-30033', 'mp-30031', 'mp-9272', 'mp-4034', 'mp-16813', 'mp-4613', 'mp-4343', 'mp-19635', 'mp-19637', 'mp-22971', 'mp-12735', 'mp-5745', 'mp-27982', 'mp-27980', 'mp-27989', 'mp-31984', 'mp-27296', 'mp-27297', 'mp-30106', 'mp-8253', 'mp-960', 'mp-642803', 'mp-22064', 'mp-23514', 'mp-3342', 'mp-12120', 'mp-29291', 'mp-21325', 'mp-29109', 'mp-25469', 'mp-29102', 'mp-31077', 'mp-31076', 'mp-31070', 'mp-20422', 'mp-27431', 'mp-27436', 'mp-27439', 'mp-21118', 'mp-28741', 'mp-28745', 'mp-7715', 'mp-20518', 'mp-20519', 'mp-604915', 'mp-675367', 'mp-20227', 'mp-27502', 'mp-740', 'mp-13051', 'mp-9479', 'mp-13292', 'mp-28256', 'mp-28253', 'mp-28251', 'mp-8652', 'mp-27616', 'mp-19912', 'mp-8890', 'mp-6598', 'mp-28414', 'mp-2928', 'mp-23192', 'mp-6913', 'mp-6912', 'mp-30938', 'mp-30937', 'mp-368', 'mp-24422', 'mp-30930', 'mp-566', 'mp-1690', 'mp-16755', 'mp-1100', 'mp-1103', 'mp-1852', 'mp-2815', 'mp-572465', 'mp-23041', 'mp-23046', 'mp-10188', 'mp-13870', 'mp-29928', 'mp-29929', 'mp-27358', 'mp-27351', 'mp-11342', 'mp-27199', 'mp-8315', 'mp-30022', 'mp-30025', 'mp-19139', 'mp-27724', 'mp-19134', 'mp-17407', 'mp-30200', 'mp-22101', 'mp-23457', 'mp-27994', 'mp-23292', 'mp-23293', 'mp-550070', 'mp-3977', 'mp-23294', 'mp-12581', 'mp-30117', 'mp-1419', 'mp-3404', 'mp-4978', 'mp-27979', 'mp-27978', 'mp-27976', 'mp-27975', 'mp-27974', 'mp-4571', 'mp-22072', 'mp-4574', 'mp-8132', 'mp-3199', 'mp-17978', 'mp-29287', 'mp-16115', 'mp-11620', 'mp-23583', 'mp-29117', 'mp-29116', 'mp-29114', 'mp-29732', 'mp-20343', 'mp-28406', 'mp-28889', 'mp-571377', 'mp-7277', 'mp-28885', 'mp-29607', 'mp-20509', 'mp-20505', 'mp-27519', 'mp-29428', 'mp-29423', 'mp-29422', 'mp-27512', 'mp-29427', 'mp-604908', 'mp-772', 'mp-7472', 'mp-10998', 'mp-10990', 'mp-28647', 'mp-11801', 'mp-28064', 'mp-28066', 'mp-28067', 'mp-7988', 'mp-6050', 'mp-29419', 'mp-18344', 'mp-28460', 'mp-28393', 'mp-28390', 'mp-28391', 'mp-12910', 'mp-14474', 'mp-6907', 'mp-2793', 'mp-396', 'mp-21916', 'mp-30908', 'mp-27771', 'mp-27770', 'mp-27772', 'mp-16762', 'mp-16763', 'mp-19838', 'mp-16766', 'mp-20321', 'mp-30107', 'mp-28268', 'mp-1443', 'mp-1667', 'mp-22813', 'mp-8415', 'mp-11028', 'mp-10213', 'mp-23079', 'mp-23072', 'mp-23071', 'mp-7461', 'mp-27362', 'mp-27361', 'mp-27367', 'mp-20885', 'mp-9254', 'mp-9251', 'mp-9252', 'mp-28092', 'mp-9095', 'mp-5589', 'mp-23440', 'mp-23288', 'mp-12597', 'mp-31473', 'mp-23281', 'mp-23280', 'mp-23282', 'mp-29709', 'mp-10030', 'mp-10033', 'mp-24945', 'mp-24944', 'mp-27948', 'mp-27949', 'mp-30302', 'mp-27943', 'mp-7794', 'mp-22534', 'mp-8124', 'mp-8129', 'mp-27386', 'mp-2798', 'mp-557719', 'mp-29123', 'mp-29124', 'mvc-13245', 'mp-22855', 'mp-21358', 'mp-7260', 'mp-28722', 'mp-28721', 'mp-22623', 'mp-23712', 'mp-22628', 'mp-22629', 'mp-8781', 'mp-7868', 'mp-762', 'mp-23844', 'mp-20713', 'mp-20242', 'mp-28095', 'mp-28094', 'mp-7462', 'mp-28096', 'mp-7466', 'mp-13274', 'mp-15800', 'mp-12857', 'mp-28099', 'mp-31000', 'mp-29708', 'mp-28650', 'mp-11815', 'mp-29254', 'mp-29257', 'mp-28070', 'mp-28077', 'mp-28074', 'mp-19933', 'mp-572892', 'mp-7991', 'mp-5824', 'mp-6023', 'mp-17172', 'mp-28384', 'mp-28474', 'mp-28470', 'mp-28479', 'mp-2278', 'mp-6930', 'mp-27783', 'mp-2273', 'mp-30911', 'mp-30913', 'mp-383', 'mp-27702', 'mp-27708', 'mp-13093', 'mp-17088', 'mp-1879', 'mp-22399', 'mp-5350', 'mp-22392', 'mp-2156', 'mp-23065', 'mp-23066', 'mp-680308', 'mp-23068', 'mp-9798', 'mp-9797', 'mp-19795', 'mp-19790', 'mp-23317', 'mp-20890', 'mp-27378', 'mp-11324', 'mp-9247', 'mp-27178', 'mp-23150', 'mp-5597', 'mp-5770', 'mp-2275', 'mp-23472', 'mp-23273', 'mp-11251', 'mp-16911', 'mp-16913', 'mp-10027', 'mp-9330', 'mp-27952', 'mp-31220', 'mp-3205', 'mp-20356', 'mp-12672', 'mp-20879', 'mp-27393', 'mp-27391', 'mp-27397', 'mp-27394', 'mp-4425', 'mp-32535', 'mp-11646', 'mp-32536', 'mp-12119', 'mp-11641', 'mp-22701', 'mp-29133', 'mp-29132', 'mp-23652', 'mp-23650', 'mp-9561', 'mp-23498', 'mp-23496', 'mp-9889', 'mp-9881', 'mp-10681', 'mp-30332', 'mp-21439', 'mp-21380', 'mp-28860', 'mp-28711', 'mp-28866', 'mp-20526', 'mp-23724', 'mp-7859', 'mvc-15971', 'mp-29403', 'mp-29402', 'mp-29407', 'mp-29408', 'mp-7853', 'mp-20253', 'mp-20091', 'mp-7459', 'mp-28667', 'mp-29249', 'mp-11860', 'mp-29718', 'mp-29713', 'mp-28049', 'mp-28046', 'mp-19921', 'mp-28448', 'mp-28449', 'mp-28442', 'mp-7038', 'mp-2207', 'mp-24416', 'mp-8752', 'mp-9650', 'mp-1468', 'mp-1469', 'mp-28282', 'mp-28283', 'mp-2668', 'mp-18785', 'mp-7149', 'mp-2841', 'mp-1319', 'mp-1863', 'mp-1644', 'mp-1861', 'mp-1641', 'mp-8436', 'mp-604910', 'mp-10232', 'mp-1797', 'mp-449', 'mp-17789', 'mp-6920', 'mp-19074', 'mp-12681', 'mp-12682', 'mp-4160', 'mp-19079', 'mp-4893', 'mp-23909', 'mp-27507', 'mp-23121', 'mp-10204', 'mp-19496', 'mp-17438', 'mp-30099', 'mp-30092', 'mp-30096', 'mp-4093', 'mp-23267', 'mp-541582', 'mp-23264', 'mp-11262', 'mp-4762', 'mp-48', 'mp-27924', 'mp-27921', 'mp-27922', 'mp-27928', 'mp-27929', 'mp-12062', 'mp-13792', 'mp-23397', 'mp-23396', 'mp-27815', 'mp-31304', 'mp-30527', 'mp-31900', 'mp-540625', 'mp-17945', 'mp-4437', 'mp-22737', 'mp-11493', 'mp-17891', 'mp-3010', 'mp-22412', 'mp-7208', 'mp-28701', 'mp-29659', 'mp-25558', 'mp-28709', 'mp-29652', 'mp-11784', 'mp-11787', 'mp-571122', 'mp-23739', 'mp-11788', 'mp-29097', 'mp-29095', 'mp-27543', 'mp-703', 'mp-20063', 'mp-28671', 'mp-28672', 'mp-15483', 'mp-28051', 'mp-8695', 'mp-29502', 'mp-29504', 'mp-30890', 'mp-995193', 'mp-580226', 'mp-634', 'mp-7534', 'mp-27490', 'mp-28458', 'mp-2076', 'mp-10651', 'mp-10656', 'mp-7028', 'mp-7029', 'mp-23866', 'mp-14444', 'mp-14338', 'mp-10532', 'mp-27727', 'mp-19862', 'mp-27725', 'mp-523', 'mp-525', 'mp-18021', 'mp-2507', 'mp-15559', 'mp-29898', 'mp-1634', 'mp-10228', 'mp-23007', 'mp-23003', 'mp-23002', 'mp-6953', 'mp-567441', 'mp-27153', 'mp-11300', 'mp-8590', 'mp-30063', 'mp-30089', 'mp-557592', 'mp-13658', 'mp-15776', 'mp-23256', 'mp-23250', 'mp-29870', 'mp-31428', 'mp-10006', 'mp-10009', 'mp-4937', 'mp-30409', 'mp-27935', 'mp-27934', 'mp-12091', 'mp-31337', 'mp-30248', 'mp-3530', 'mp-1213', 'mp-3312', 'mp-1214', 'mp-4449', 'mp-17939', 'mp-4445', 'mp-572758', 'mp-31406', 'mp-21413', 'mp-21410', 'mp-3006', 'mp-28693', 'mp-7237', 'mp-12264', 'mp-28846', 'mp-10144', 'mp-28845', 'mp-21096', 'mp-29081', 'mp-29467', 'mp-10798', 'mp-29465', 'mp-27553', 'mp-29469', 'mp-8911', 'mp-9429', 'mp-28608', 'mp-28606', 'mp-28601', 'mp-29263', 'mp-29264', 'mp-20476', 'mp-20470', 'mp-29519', 'mp-23925', 'mp-29510', 'mp-2447', 'mp-2194', 'mp-20617', 'mp-6789', 'mp-3208', 'mp-602', 'mp-6781', 'mp-605', 'mp-17146', 'mp-27485', 'mp-10648', 'mp-13134', 'mp-18305', 'mp-10643', 'mp-570062', 'mp-19985', 'mp-10645', 'mp-6543', 'mp-7784', 'mp-28116', 'mp-28117', 'mp-28118', 'mp-7839', 'mp-28330', 'mp-9631', 'mp-27739', 'mp-19878', 'mp-27734', 'mp-20107', 'mp-20106', 'mvc-13985', 'mp-20104', 'mp-1404', 'mp-18010', 'mp-29881', 'mp-29882', 'mp-10259', 'mp-1754', 'mp-1752', 'mp-5096', 'mp-8585', 'mp-19278', 'mp-12794', 'mp-23244', 'mp-23247', 'mp-27230', 'mp-29863', 'mp-29862', 'mp-11209', 'mp-10074', 'mp-3496', 'mp-32497', 'mp-12085', 'mp-6462', 'mp-10101', 'mp-3813', 'mp-17924', 'mp-22753', 'mp-984', 'mp-571033', 'mp-31037', 'mp-28938', 'mp-980', 'mp-28934', 'mp-13682', 'mp-28933', 'mp-13680', 'mp-3257', 'mp-3255', 'mp-27902', 'mp-21409', 'mp-567279', 'mp-21405', 'mp-9984', 'mp-10810', 'mp-3075', 'mp-13357', 'mp-29072', 'mp-29073', 'mp-32539', 'mp-23759', 'mp-720', 'mp-29694', 'mp-726', 'mp-24290', 'mp-20597', 'mp-21062', 'mp-28583', 'mp-28580', 'mp-12324', 'mp-12325', 'mp-28616', 'mp-22574', 'mp-22577', 'mp-22573', 'mp-571210', 'mp-28037', 'mp-20629', 'mp-612', 'mp-7954', 'mp-7487', 'mp-7483', 'mp-28038', 'mp-16403', 'mp-13126', 'mp-10678', 'mp-10676', 'mp-7559', 'mp-6573', 'mp-29458', 'mp-28100', 'mp-7826', 'mp-29453', 'mp-23801', 'mp-29327', 'mp-28329', 'mp-28328', 'mp-19431', 'mp-28326', 'mp-502', 'mp-9622', 'mp-500', 'mp-19846', 'mp-20138', 'mp-19845', 'mp-20134', 'mp-19849', 'mp-20130', 'mp-8765', 'mp-2450', 'mp-1186', 'mp-15573', 'mp-29843', 'mvc-14413', 'mp-23020', 'mp-30058', 'mp-470', 'mp-2232', 'mp-2231', 'mp-19990', 'mp-2237', 'mp-22998', 'mp-30023', 'mp-1507', 'mp-1760', 'mp-29988', 'mp-1984', 'mp-8378', 'mp-9281', 'mp-8374', 'mp-6999', 'mp-18868', 'mp-23112', 'mp-8573', 'mvc-6946', 'mvc-6948', 'mp-12523', 'mp-23238', 'mp-23230', 'mp-11217', 'mp-20902', 'mp-9371', 'mp-11213', 'mp-30427', 'mp-24936', 'mp-29853', 'mp-10065', 'mp-8482', 'mp-4138', 'mp-5625', 'mp-14815', 'mp-8677', 'mp-24809', 'mp-541449', 'mp-21273', 'mp-3667', 'mp-4663', 'mp-4662', 'mp-12748', 'mp-23697', 'mp-23694', 'mp-28929', 'mp-28928', 'mp-27910', 'mp-29196', 'mp-29194', 'mp-29190', 'mp-21474', 'mp-9996', 'mp-29666', 'mp-29663', 'mp-29067', 'mp-30999', 'mp-23563', 'mp-29689', 'mp-8976', 'mp-8140', 'mp-21055', 'mp-23355', 'mp-21188', 'mp-20459', 'mp-21183', 'mp-572822', 'mp-665', 'mvc-14734', 'mp-7962', 'mp-7964', 'mp-23297', 'mp-8806', 'mp-30952', 'mp-28483', 'mp-28480', 'mp-28136', 'mp-28135', 'mp-28130', 'mp-28131', 'mp-7814', 'mp-29337', 'mp-29440', 'mp-29443', 'mp-7760', 'mp-349', 'mp-7582', 'mp-28319', 'mp-28313', 'mp-20121', 'mp-20120', 'mp-19857', 'mp-20122', 'mp-20125', 'mp-20128', 'mp-6368', 'mp-18075', 'mp-20315', 'mp-22996', 'mp-28553', 'mp-20830', 'mp-25037', 'mp-13976', 'mp-2556', 'mp-7617', 'mp-22987', 'mp-19984', 'mp-22983', 'mp-7618', 'mp-12406', 'mp-22969', 'mp-22963', 'mp-22966', 'mp-22965', 'mp-23164', 'mp-30988', 'mp-23162', 'mp-19187', 'mp-23169', 'mp-30984', 'mp-1602', 'mp-5328', 'mp-1605', 'mp-19216', 'mp-15742', 'mp-19433', 'mp-23229', 'mp-23228', 'mp-23226', 'mp-23225', 'mp-23223', 'mp-10058', 'mp-27215', 'mp-27213', 'mp-571442', 'mp-9386', 'mp-9385', 'mp-9384', 'mp-27219', 'mp-7558', 'mp-9146', 'mp-29361', 'mp-19727', 'mp-23353', 'mp-3871', 'mp-27859', 'mp-12621', 'mp-27854', 'mp-570880', 'mp-27857', 'mp-27850', 'mp-24810', 'mp-30295', 'mp-24814', 'mp-10128', 'mp-30299', 'mp-4674', 'mp-22777', 'mp-23680', 'mp-4678', 'mp-28361', 'mp-12751', 'mp-12750', 'mp-22863', 'mp-30362', 'mp-21468', 'mp-31015', 'mp-567809', 'mp-28838', 'mp-28836', 'mp-28368', 'mp-28832', 'mvc-13391', 'mp-16053', 'mp-32450', 'mp-29057', 'mp-21260', 'mp-21264', 'mp-21267', 'mp-20026', 'mp-21688', 'mp-10900', 'mp-32131', 'mp-3123', 'mp-24030', 'mp-10906', 'mp-22558', 'mp-28633', 'mp-29236', 'mp-29234', 'mp-571279', 'mp-27697', 'mp-29543', 'mp-27699', 'mp-20333', 'mp-27455', 'mp-9810', 'mp-28491', 'mp-28490', 'mp-9548', 'mp-20480', 'mp-10614', 'mp-28128', 'mp-10616', 'mp-22601', 'mp-28121', 'mp-29308', 'mp-29309', 'mp-7775', 'mp-22424', 'mp-29305', 'mp-29300', 'mp-28301', 'mp-7595', 'mp-28304', 'mp-28306', 'mp-28308', 'mp-9600', 'mp-10728', 'mp-29568', 'mp-6353', 'mp-28011', 'mp-28013', 'mp-28544', 'mp-13963', 'mp-1168', 'mp-9776', 'mp-675438', 'mp-1563', 'mp-1566', 'mp-1568', 'mp-7609', 'mp-14112', 'mp-22959', 'mp-23177', 'mp-23176', 'mp-306', 'mp-23173', 'mp-11181', 'mp-23170', 'mp-30993', 'mp-23178', 'mp-22889', 'mp-19888', 'mp-22880', 'mp-22881', 'mp-22883', 'mp-5045', 'mp-23210', 'mp-23216', 'mp-23219', 'mp-11232', 'mp-22310', 'mp-9399', 'mp-9390', 'mp-9391', 'mp-31484', 'mp-31487', 'mp-3468', 'mp-32531', 'mp-31489', 'mp-29866', 'mp-19715', 'mp-5603', 'mp-1943', 'mp-11524', 'mp-27849', 'mp-27848', 'mp-29940', 'mp-1948', 'mp-30282', 'mp-30287', 'mp-30286', 'mp-24866', 'mp-4488', 'mp-2129', 'mp-568592', 'mp-4649', 'mp-3592', 'mp-30009', 'mp-28947', 'mp-22161', 'mp-28944', 'mp-12769', 'mp-27656', 'mp-4246', 'mp-16789', 'mp-4249', 'mp-8045', 'mp-3047', 'mp-685055', 'mp-10863', 'mp-31020', 'mp-32491', 'mp-30341', 'mp-29048', 'mp-23546', 'mp-21210', 'mp-542096', 'mp-29238', 'mp-8953', 'mp-10913', 'mp-21037', 'mp-21038', 'mp-12351', 'mp-11887', 'mp-29555', 'mp-27683', 'mp-27682', 'mp-1115', 'mp-641', 'mp-21149', 'mp-27441', 'mp-27440', 'mp-13178', 'mp-20494', 'mp-20498', 'mp-27449', 'mp-8825', 'mp-13173', 'mp-12992', 'mp-28153', 'mp-29311', 'mp-29310', 'mp-29314', 'mp-7293', 'mp-28791', 'mp-7744', 'mp-28374', 'mp-569715', 'mp-28372', 'mp-14131', 'mp-8732', 'mp-794', 'mp-7615', 'mp-1683', 'mp-28577', 'mp-28575', 'mp-28572', 'mp-28570', 'mp-28208', 'mp-1170', 'mp-29541', 'mp-1572', 'mp-27662', 'mp-27666', 'mp-27664', 'mp-7634', 'mp-8435', 'mp-29420', 'mp-12422', 'mp-22945', 'mp-2977', 'mp-23141', 'mp-10231', 'mp-9010', 'mp-1379', 'mp-22897', 'mp-22893', 'mp-22892', 'mp-5073', 'mp-23205', 'mp-23204', 'mp-23207', 'mp-1397', 'mp-27277', 'mp-19524', 'mp-19527', 'mp-31495', 'mp-10615', 'mp-13774', 'mp-23375', 'mp-23376', 'mp-27302', 'mp-27871', 'mp-11553', 'mp-29950', 'mp-24877', 'mp-10145', 'mp-2136', 'mp-541183', 'mp-11433', 'mp-5709', 'mp-12770', 'mp-12776', 'mp-557500', 'mp-3922', 'mp-3700', 'mp-21335', 'mp-31036', 'mp-31038', 'mp-21339', 'mp-9948', 'mp-11722', 'mp-849079', 'mp-23553', 'mp-21208', 'mp-21204', 'mp-20551', 'mp-20004', 'mp-21020', 'mp-32541', 'mp-9830', 'mp-27475', 'mp-18609', 'mp-20311', 'mp-29496', 'mp-22446', 'mp-20797', 'mp-29491', 'mp-7280', 'mp-22449', 'mp-29498', 'mp-27788', 'mp-12948', 'mp-27785', 'mp-27780', 'mp-20177', 'mp-28369', 'mp-20175', 'mp-20779', 'mp-17004', 'mp-18731', 'mp-2418', 'mp-28569', 'mp-2145', 'mp-2416', 'mp-29567', 'mp-2414', 'mp-23970', 'mp-8708', 'mp-1492', 'mp-28211', 'mp-28214', 'mp-430', 'mp-27659', 'mp-8616', 'mp-27655', 'mp-27652', 'mp-8613', 'mp-7622', 'mp-5282', 'mp-505366', 'mp-34289', 'mp-22932', 'mp-22939', 'mp-23386', 'mp-30974', 'mp-24468', 'mp-30979', 'mp-30978', 'mp-2565', 'mp-6669', 'mp-25470', 'mp-22849', 'mp-14623', 'mp-29671', 'mp-224', 'mp-1814', 'mp-29817', 'mp-18958', 'mp-30770', 'mp-30777', 'mp-5443', 'mp-5446', 'mp-23309', 'mp-23307', 'mp-24848', 'mp-27866', 'mp-27863', 'mp-8190', 'mp-27319', 'mp-27868', 'mp-4004', 'mp-9228', 'mp-13498', 'mp-13497', 'mp-568846', 'mp-22144', 'mp-28963', 'mp-23412', 'mp-12707', 'mp-8215', 'mp-31045', 'mp-31040', 'mp-9950', 'mp-13383', 'mp-13385', 'mp-4530', 'mp-13542', 'mp-29028', 'mp-29027', 'mp-13548', 'mp-29020', 'mp-23521', 'mp-21239', 'mp-9481', 'mp-12376', 'mp-12372', 'mp-29797', 'mp-29154', 'mp-8845', 'mp-28683', 'mp-27461', 'mp-8848', 'mp-27465', 'mp-20307', 'mp-29774', 'mp-13152', 'mp-10628', 'mp-9828', 'mp-9827', 'mp-29489', 'mp-28176', 'mp-28174', 'mp-28175', 'mp-28178', 'mp-1118', 'mp-11973', 'mp-20169', 'mp-28192', 'mp-20762', 'mp-20765', 'mp-20764', 'mp-24204', 'mp-24205', 'mp-13024', 'mp-8999', 'mp-8998', 'mp-7925', 'mp-7926', 'mp-25078', 'mp-29573', 'mp-2355', 'mp-23099', 'mp-28228', 'mp-6384', 'mp-28223', 'mp-28220', 'mp-27641', 'mp-19941', 'mp-27645', 'mp-27647', 'mp-27648', 'mp-6708', 'mp-17708', 'mp-28429', 'mp-28358', 'mp-6096', 'mp-6094', 'mp-28421', 'mp-28357', 'mp-28426', 'mp-22923', 'mp-30946', 'mp-358', 'mp-30945', 'mp-30943', 'mp-17637', 'mp-617', 'mp-5014', 'mp-27253', 'mp-22859', 'mp-235', 'mp-27257', 'mp-236', 'mp-22853', 'mp-22850', 'mp-23194', 'mp-14659', 'mp-29803', 'mp-22854', 'mp-9102', 'mp-562100', 'mp-29526', 'mp-5470', 'mp-17524', 'mp-30125', 'mp-4391', 'mp-19761', 'mp-23318', 'mp-15679', 'mp-23316', 'mp-23315', 'mp-23313', 'mp-23312', 'mp-22210', 'mp-4384', 'mp-22217', 'mp-6473', 'mp-30018', 'mp-3891', 'mp-20747', 'mp-28996', 'mp-31725', 'mp-541070', 'mp-20994', 'mp-9922', 'mp-9921', 'mp-9920', 'mp-31053', 'mp-909', 'mp-17835', 'mp-11701', 'mp-29016', 'mp-29015', 'mp-29018', 'mp-23536', 'mp-3387', 'mp-7008', 'mp-3382', 'mp-567832', 'mp-1016', 'mp-13384', 'mp-20579', 'mp-13449', 'mp-29459', 'mp-27418', 'mp-696944', 'mp-28691', 'mp-28919', 'mp-29744', 'mp-27411', 'mp-18625', 'mp-13146', 'mp-29349', 'mp-570219', 'mp-971787', 'mp-10897', 'mp-28765', 'mp-29026', 'mp-22661', 'mp-8390', 'mp-22667', 'mp-28189', 'mp-20195', 'mp-20198', 'mp-625199', 'mp-20203', 'mp-28916', 'mp-27529', 'mp-20754', 'mp-10761', 'mp-2160', 'mp-28500', 'mp-7932', 'mp-7425', 'mp-28509', 'mp-28230', 'mp-28233', 'mp-5864', 'mp-10426', 'mp-9731', 'mp-10421', 'mp-29831', 'mp-27634', 'mp-27630', 'mp-18156', 'mp-18407', 'mp-694', 'mp-27639', 'mp-690', 'mp-28341', 'mp-7645', 'mp-1581', 'mp-19372', 'mp-1582', 'mp-567354', 'mp-22917', 'mp-30953', 'mp-8516', 'mp-30954', 'mp-10553', 'mp-341', 'mp-27743', 'mp-10326', 'mp-32887', 'mp-569581', 'mp-1125', 'mp-12546', 'mp-27245', 'mp-27246', 'mp-27240', 'mp-22868', 'mp-22869', 'mp-27249', 'mp-1834', 'mp-22861', 'mp-30113', 'mp-19578', 'mp-4731', 'mp-30139', 'mp-19755', 'mp-10838', 'mp-23324', 'mp-5489', 'mp-27337', 'mp-10177', 'mp-11364', 'mp-30003', 'mp-10371', 'mp-30228', 'mp-3886', 'mp-30224', 'mp-28983', 'mp-28982', 'mp-22125', 'mp-23434', 'mp-23435', 'mp-25667', 'mp-3421', 'mp-31999', 'mp-31268', 'mp-3779', 'mp-28878', 'mp-9936', 'mp-974', 'mp-22691', 'mp-22693', 'mp-23504', 'mp-3392', 'mp-31198', 'mp-21563', 'mp-804', 'mp-11786', 'mp-17287', 'mp-22582', 'mp-29179', 'mp-29174', 'mp-29941', 'mp-21363', 'mp-9592', 'mp-31474', 'mp-21106', 'mp-21103', 'mp-27399', 'mp-7708', 'mp-7705', 'mp-7701', 'mp-28757', 'mp-28754', 'mp-28753', 'mp-20560', 'mp-27462', 'mp-20184', 'mp-20210', 'mp-20215', 'mp-29593', 'mp-29592', 'mp-29599', 'mp-11826', 'mp-18637', 'mp-28247', 'mp-569623', 'mp-28248', 'mp-2336', 'mp-8644', 'mp-22992', 'mp-18169', 'mp-27628', 'mp-7077', 'mp-972794', 'mp-2755', 'mp-28407', 'mp-28405', 'mp-6070', 'mp-27944', 'mp-37514', 'mp-22908', 'mp-1224', 'mp-505531', 'mp-23184', 'mp-23185', 'mp-23189', 'mp-7724', 'mp-10569', 'mp-13995', 'mp-1136', 'mp-6615']
import sys
print len(new),'len'
#sys.exit()

def make_subscript_spg_html(string=''):
    orig=string
    for i,char in enumerate(string):
        if char=='-':
           char=string[i+1]
           new_char=str(char)+str('&#773;')#str("$")+str("\\")+str("bar")+str("{")+str(char)+str("}$")#r'$\alpha_i
           #new_char=str(char)+str('&#772;')#str("$")+str("\\")+str("bar")+str("{")+str(char)+str("}$")#r'$\alpha_i
           #new_char=str(char)+str('&oline;')#str("$")+str("\\")+str("bar")+str("{")+str(char)+str("}$")#r'$\alpha_i
           #new_char=str("'$")+str("_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str("r'$")+str("\_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str('_')+str(char)
           string=string.replace(char,new_char)
           string=string.replace('-','')
        if char=='_' and 'sub' not in string:
           char=string[i+1]
           new_char=str("<sub>")+str(string[i+1])+str("</sub>")#r'$\alpha_i
           #new_char=str("<sub>")+str(char)+str("</sub>")#r'$\alpha_i
           #new_char=str("'$")+str("_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str("r'$")+str("\_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str('_')+str(char)
           string=string.replace(char,new_char)
           string=string.replace('_','')
           #print orig,string
    return string

def make_subscript_redf_html(string=''):
    for i,char in enumerate(string):
        if char.isdigit() :
           new_char=str("<sub>")+str(char)+str("</sub>")#r'$\alpha_i
           #new_char=str("'$")+str("_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str("r'$")+str("\_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str('_')+str(char)
           string=string.replace(char,new_char)
    return string



def make_subscript_spg(string=''):
    for i,char in enumerate(string):
        if char=='-':
           char=string[i]+string[i+1]
           new_char=str("$")+str("\\")+str("bar")+str("{")+str(char)+str("}$")#r'$\alpha_i
           #new_char=str("'$")+str("_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str("r'$")+str("\_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str('_')+str(char)
           string=string.replace(char,new_char)
           string=string.replace('-','')
        if char=='_':
           char=string[i]+string[i+1]
           new_char=str("$")+str(char)+str("$")#r'$\alpha_i
           #new_char=str("'$")+str("_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str("r'$")+str("\_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str('_')+str(char)
           string=string.replace(char,new_char)
           string=string.replace('-','')
    return string

def make_subscript_redf(string=''):
    for i,char in enumerate(string):
        if char.isdigit() :
           new_char=str("$")+str("_")+str(char)+str("$")#r'$\alpha_i
           #new_char=str("'$")+str("_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str("r'$")+str("\_")+str(char)+str("$'")#r'$\alpha_i
           #new_char=str('_')+str(char)
           string=string.replace(char,new_char)
    string=string.replace('$_$','$')
    string=string.replace('$$','$')
    return string

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.1f}%  '.format(p=pct)
        #return '{p:.1f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def make_pie_dict(data=None,name='pie.png',select=7,grid=None):
        keys=[]
        values=[]
        for k,v in data.iteritems():
           keys.append(k)
           values.append(v)
        order = np.argsort(values)
        order=order[::-1]
        
        keys=np.array(keys)[order]
        values=np.array(values)[order]

        new_keys=keys[0:select]
        new_values=values[0:select]
        colors=['yellowgreen', 'gold', 'lightskyblue', 'lightcoral','lightblue','lightgreen','olive','tomato']
        a=new_keys
        b=new_values
        if 'Crys' not in name:
           a=np.append(new_keys,'others')
           sum=np.sum(values[select:len(values)])
           b=np.append(new_values,sum)
        c=[]
        if 'Ano' in name:
          for i in a:
            c.append(make_subscript_redf(i))
          print 'c for aAno',c
          #c=a
        if 'Spg' in name:
          for i in a:
            c.append(make_subscript_spg(i))
        if len(c)!=len(a):
            c=a
        else:
            print 'len c a',len(c),len(a)
        #b=np.append(new_values,np.sum(keys[select:len(keys)]))
        #print a,b,c
        plt.subplot(grid, aspect=1)
        #print keys,values
        #plt = get_publication_quality_plot(12, 14)
        plt.pie(b,labels=c, autopct=make_autopct(b),pctdistance=0.8,colors=colors)
        #ax1.plot(eg,n, 'o')
        #ax1.set_aspect(1)
        #ax1.plot(mp,mp)
        #plt.tight_layout()
        #filename=name
        #plt.savefig(filename)
        #plt.close()
def main_func(perc=5,new=[],data=[]):
       count=0
       #f=open('info.json','r')
       #import json
       #data=json.load(f,cls=MontyDecoder)
       #f.close()

       mem={}
       memspg={}
       memnary={}
       memcrys={}
       memabc={}
       mem1={}
       memspg1={}
       memnary1={}
       memcrys1={}
       memabc1={}
       chalcogen=[] #O,S,Se,Te,Po
       halogen=[]#F,Cl,Br,I,At
       pnictogen=[]#N,P,As,Sb,Bi
       othergen=[]
       aa=[]
       bb=[]
       cc=[]
       for i in data:
         mmpid= str(i['mp_id'])
         if mmpid in new:

               ini= (i['ini_structure'])#Poscar(i['ini'][0]).structure
               fin= (i['fin_structure'])#Poscar(i['fin'][0]).structure
               a,b,c=ini.lattice.abc
               a1,b1,c1=fin.lattice.abc
               ratioc= 100*round(abs((float(c)-float(c1))/float(c1)),2)
               ratiob= 100*round(abs((float(b)-float(b1))/float(b1)),2)
               ratioa= 100*round(abs((float(a)-float(a1))/float(a1)),2)
               sp=len(ini.types_of_specie)
               icsd=i["icsd"]
               mula=str(ini.composition.reduced_formula)
               rat='na'
               if  not (round(a1,3)==round(b1,3)==round(c1,3))  and sp>1 and icsd!=None and icsd!=[] :
                  if  ratioc>=perc :
                       
                      rat='c'
                      aa.append(mmpid)
                  if  ratiob>=perc :
                      rat='b'
                      bb.append(mmpid)
                  if  ratioa>=perc :
                      rat='a'
                      cc.append(mmpid)
                    
               count=count+1
               count=count+1
               crys=str(SpacegroupAnalyzer(ini).get_crystal_system())
               sgp=str(SpacegroupAnalyzer(ini).get_spacegroup_symbol())
               anonymized_formula=str(ini.composition.anonymized_formula)
               if anonymized_formula not in mem:
                      mem[anonymized_formula]=1
               else:
                      mem[anonymized_formula]=mem[anonymized_formula]+1
               if anonymized_formula not in mem1:
                      mem1[anonymized_formula]=mula
               else:
                      mem1[anonymized_formula]=str(mem1[anonymized_formula])+str(', ')+str(mula)



               if sgp not in memspg:
                      memspg[sgp]=1
               else:
                      memspg[sgp]=memspg[sgp]+1
               if sgp not in memspg1:
                      memspg1[sgp]=mula
               else:
                      memspg1[sgp]=str(memspg1[sgp])+str(', ')+str(mula)



               if rat not in memabc:
                      memabc[rat]=1
               else:
                      memabc[rat]=memabc[rat]+1
               if rat not in memabc1:
                      memabc1[rat]=1
               else:
                      memabc1[rat]=memabc1[rat]+1




               if crys not in memcrys:
                      memcrys[crys]=1
               else:
                      memcrys[crys]=memcrys[crys]+1
               if crys not in memcrys1:
                      memcrys1[crys]=mula
               else:
                      memcrys1[crys]=str(memcrys1[crys])+str(', ')+str(mula)




               nary=''
               if len(ini.composition.elements)==1:
                   nary='Unary'
               elif len(ini.composition.elements)==2:
                   nary='Binary'
               elif len(ini.composition.elements)==3:
                   nary='Ternary'
               elif len(ini.composition.elements)==4:
                   nary='Quaternary'
               elif len(ini.composition.elements)==5:
                   nary='Quinary'
               elif len(ini.composition.elements)==6:
                   nary='Senary'
               elif len(ini.composition.elements)==7:
                   nary='Septenary'
               else:
                   nary='MultiComponent'


               if nary not in memnary:
                      memnary[nary]=1
               else:
                      memnary[nary]=memnary[nary]+1
               p=ini

               if any(i in p.composition.elements for i in [Element('O'),Element('S'),Element('Se'),Element('Te'),Element('Po')]):
                  chalcogen.append(mmpid)
               if any(i in p.composition.elements for i in [Element('F'),Element('Cl'),Element('Br'),Element('I'),Element('At')]):
                  halogen.append(mmpid)
               if any(i in p.composition.elements for i in [Element('N'),Element('P'),Element('As'),Element('Sb'),Element('Bi')]):
                   pnictogen.append(mmpid)

               if not any(i in p.composition.elements for i in [Element('N'),Element('P'),Element('As'),Element('Sb'),Element('Bi'),Element('F'),Element('Cl'),Element('Br'),Element('I'),Element('At'),Element('O'),Element('S'),Element('Se'),Element('Te'),Element('Po')]):
                  othergen.append(mmpid)


       plt = get_publication_quality_plot(12, 14)
       set1 = set(chalcogen)
       set2 = set(halogen)
       set3 = set(pnictogen)
       seta = set(aa)
       setb = set(bb)
       setc = set(cc)
       #print set1
       #print set2
       #print set3,len(data)

       








       plt.subplot(the_grid[0, 0], aspect=1)
       venn3([seta, setb, setc], ('latt-a', 'latt-b', 'latt-c'))
       plt.title('(a)')
       plt.subplot(the_grid[0, 1], aspect=1)
       venn3([set1, set2, set3], ('Chalcogenide', 'Halide', 'Pnictide'))
       plt.title('(b)')
       make_pie_dict(data=mem,name='AnoPie.png',grid=the_grid[1,0])
       plt.title('(c)')
       make_pie_dict(data=memspg,name='SpgPie.png',grid=the_grid[1,1])
       plt.title('(d)')
       make_pie_dict(data=memcrys,name='CrysPie.png',select=6,grid=the_grid[2,0])
       plt.title('(e)')
       make_pie_dict(data=memnary,name='NaryPie.png',select=2,grid=the_grid[2,1])
       plt.title('(f)')
       name=str('Venn-')+str(perc)+str('.png')
       plt.savefig(name)
       plt.close()
       print "count",count,len(new)#for k,v in data.iteritems():
       if perc==5:
		h=open('Ano.html','w')
		line=str("      <style>")+'\n'
		h.write(line)
		line=str('table, td, th {')+'\n'
		h.write(line)
		line=str('    border: 1px solid black;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)

		line=str('table {')+'\n'
		h.write(line)
		line=str('    border-collapse: collapse;')+'\n'
		h.write(line)
		line=str('    width: 30%;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)

		line=str('th {')+'\n'
		h.write(line)
		line=str('    height: 50px;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)
		line=str('</style>')+'\n'
		h.write(line)

		line=str('<table >')+'\n'
                h.write(line)
                for i,j in mem1.iteritems():
                  line=str("<tr>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(make_subscript_redf_html(i))+str("</td>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(j)+str("</td>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(len(str(j).split(',')))+str("</td>")+'\n'
                  h.write(line)
                  line=str("</tr>")+'\n'
                  h.write(line)
                line=str("</table>")+'\n'
                h.write(line)
                h.close()



#def make_subscript_spg_html(string=''):
#def make_subscript_redf_html(string=''):
		h=open('Spg.html','w')
		line=str("      <style>")+'\n'
		h.write(line)
		line=str('table, td, th {')+'\n'
		h.write(line)
		line=str('    border: 1px solid black;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)

		line=str('table {')+'\n'
		h.write(line)
		line=str('    border-collapse: collapse;')+'\n'
		h.write(line)
		line=str('    width: 30%;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)

		line=str('th {')+'\n'
		h.write(line)
		line=str('    height: 50px;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)
		line=str('</style>')+'\n'
		h.write(line)

		line=str('<table >')+'\n'
                h.write(line)
                for i,j in memspg1.iteritems():
                  line=str("<tr>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(make_subscript_spg_html(i))+str("</td>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(j)+str("</td>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(len(str(j).split(',')))+str("</td>")+'\n'
                  h.write(line)
                  line=str("</tr>")+'\n'
                  h.write(line)
                line=str("</table>")+'\n'
                h.write(line)

		h=open('Crys.html','w')
		line=str("      <style>")+'\n'
		h.write(line)
		line=str('table, td, th {')+'\n'
		h.write(line)
		line=str('    border: 1px solid black;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)

		line=str('table {')+'\n'
		h.write(line)
		line=str('    border-collapse: collapse;')+'\n'
		h.write(line)
		line=str('    width: 30%;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)

		line=str('th {')+'\n'
		h.write(line)
		line=str('    height: 50px;')+'\n'
		h.write(line)
		line=str('}')+'\n'
		h.write(line)
		line=str('</style>')+'\n'
		h.write(line)

		line=str('<table >')+'\n'
                h.write(line)
                for i,j in memcrys1.iteritems():
                  line=str("<tr>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(i)+str("</td>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(j)+str("</td>")+'\n'
                  h.write(line)
                  line=str("<td>")+str(len(str(j).split(',')))+str("</td>")+'\n'
                  h.write(line)
                  line=str("</tr>")+'\n'
                  h.write(line)
                line=str("</table>")+'\n'
                h.write(line)

       if perc==5:
         for i,j in mem1.iteritems():
            print i,j
         for i,j in memspg1.iteritems():
            print i,j
         for i,j in memcrys1.iteritems():
            print i,j


data = loadfn('/home/knc6/bin/MPall_datacopy.json', cls=MontyDecoder)
new3=[]
new5=[]
new7=[]
new10=[]
#3116 1356 819 375
for i in data:
           mmpid= str(i['mp_id'])
           ini= (i['ini_structure'])#Poscar(i['ini'][0]).structure
           fin= (i['fin_structure'])#Poscar(i['fin'][0]).structure
           a,b,c=ini.lattice.abc
           a1,b1,c1=fin.lattice.abc
           ratioc= 100*round(abs((float(c)-float(c1))/float(c1)),2)
           ratiob= 100*round(abs((float(b)-float(b1))/float(b1)),2)
           ratioa= 100*round(abs((float(a)-float(a1))/float(a1)),2)
           sp=len(ini.types_of_specie)
           icsd=i["icsd"]
           if  not (round(a1,3)==round(b1,3)==round(c1,3))  and sp>1 and icsd!=None and icsd!=[] :
             if  ratioc>=3 or ratiob>=3 or ratioa>=3:
                 new3.append(mmpid)
             if  ratioc>=5 or ratiob>=5 or ratioa>=5:
                 new5.append(mmpid)
             if  ratioc>=7 or ratiob>=7 or ratioa>=7:
                 new7.append(mmpid)
             if  ratioc>=10 or ratiob>=10 or ratioa>=10:
                 new10.append(mmpid)
print len(new3),len(new5),len(new7),len(new10)
main_func(new=new3,perc=3,data=data)
main_func(new=new5,perc=5,data=data)
main_func(new=new7,perc=7,data=data)
main_func(new=new10,perc=10,data=data)
