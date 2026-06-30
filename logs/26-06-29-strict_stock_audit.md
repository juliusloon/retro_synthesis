# Strict Buyable Stock 审计报告

生成日期: 2026-06-30 09:51

## 参数

- 源文件: `flavonoid_stock_inchikeys.txt`
- MW 截止: 350.0
- 重原子数截止: 25

## 统计

| 指标 | 数值 |
|---|---:|
| 原始条目数 | 297 |
| 保留条目数 | 261 |
| 排除条目数 | 36 |
| 未验证条目数 | 0 |

## 处理策略

- **保留**: MW ≤ 350 且重原子数 ≤ 25 的可购买小分子
- **排除**: MW 或重原子数超标的复杂天然产物/糖苷
- **未验证**: 缺少分子信息的 InChIKey，**不纳入** strict stock（单独存放于 `unverified_stock_inchikeys.txt`）

## 排除条目

| InChIKey | SMILES | MW | 重原子数 | 排除原因 |
|---|---|---:|---:|---|
| `ABYDWDZTWUDKPU-UHFFFAOYSA-N` | `CCCC[Te+](CCCC)CC(=O)c1ccccc1` | 361.0 | 18 | MW=361.0 > 350.0 |
| `ALEHFXXYPFUXSF-UHFFFAOYSA-N` | `O=S(=O)(C1=CCOc2cc(OCc3ccccc3)ccc21)c1ccccc1` | 378.4 | 27 | MW=378.4 > 350.0 |
| `AXJCTLVJVOXEIQ-BIFWDUKVSA-N` | `COc1ccc([C@@H]2CC(=O)c3c(O)cc(O[C@H]4O[C@@H](CO)[C@H](O)[C@@...` | 610.6 | 43 | MW=610.6 > 350.0 |
| `BUBVLQDEIIUIQG-LZQYHZOCSA-N` | `O=C1O[C@H](COCc2ccccc2)C(OCc2ccccc2)[C@@H](OCc2ccccc2)[C@H]1...` | 538.6 | 40 | MW=538.6 > 350.0 |
| `BWKJTESGMBODJZ-DWTLZWDMSA-N` | `CC(=O)Oc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O[C@@H]3OCC(OC(C)=O)[C...` | 1018.9 | 74 | MW=1018.9 > 350.0 |
| `CJYWNEOODUJSQE-UHFFFAOYSA-N` | `O=c1cc(COc2ccccc2I)oc2ccccc12` | 378.2 | 20 | MW=378.2 > 350.0 |
| `CYAYKKUWALRRPA-LYPMUSHJSA-N` | `CC(=O)OC[C@H]1OC(Br)[C@H](OC(C)=O)[C@H](OC(C)=O)C1OC(C)=O` | 411.2 | 24 | MW=411.2 > 350.0 |
| `CYEWBUBDSRRSRF-UHFFFAOYSA-N` | `CCOC(=O)CC(=O)c1cc([N+](=O)[O-])c([N+](=O)[O-])c([N+](=O)[O-...` | 372.2 | 26 | MW=372.2 > 350.0 |
| `GKXMWFAVHWCSJP-UHFFFAOYSA-N` | `CCCC[Sn](CCCC)(CCCC)c1ccc(C)cc1` | 381.2 | 20 | MW=381.2 > 350.0 |
| `HFDSRPIZLNFESV-TZKAHBFUSA-N` | `C[C@H]1O[C@H](Oc2cc(O)c3c(=O)c(O[C@H]4O[C@@H](CO)[C@H](O[C@H...` | 918.8 | 65 | MW=918.8 > 350.0 |
| `JBCZWJTVGBEWIR-UHFFFAOYSA-N` | `COc1cc(O)c2c(c1)OC(c1ccc(O)c(OC)c1)CC21SCCCS1` | 406.5 | 27 | MW=406.5 > 350.0 |
| `JJXATNWYELAACC-LYPMUSHJSA-N` | `CC(=O)OC[C@H]1OC(F)[C@H](OC(C)=O)[C@H](OC(C)=O)C1OC(C)=O` | 350.3 | 24 | MW=350.3 > 350.0 |
| `JPSHKIIZDYVUNB-XBUOOGCBSA-O` | `COC1=CC(=CC(CO)Oc2ccc([C@H]3Oc4cc(O)cc(O)c4C(=O)[C@@H]3O)cc2...` | 483.4 | 35 | MW=483.4 > 350.0 |
| `JSEBZQXYIRIWMB-NNZPWHBESA-N` | `CC(=O)OC[C@H]1CC(C(=O)NC(Cl)(Cl)Cl)[C@H](OC(C)=O)[C@H](OC(C)...` | 490.7 | 30 | MW=490.7 > 350.0 |
| `KQCMATLOONBARE-UHFFFAOYSA-N` | `O=C(COc1ccccc1-c1coc2ccccc2c1=O)On1ccccc1=S` | 405.4 | 29 | MW=405.4 > 350.0 |
| `LCMPHIATIKTQMI-UHFFFAOYSA-K` | `CC(=O)O[Pb](OC(C)=O)(OC(C)=O)c1ccccc1` | 461.4 | 19 | MW=461.4 > 350.0 |
| `LIILCIIZVMZTSC-UHFFFAOYSA-N` | `CCCC[Sn](CCCC)(CCCC)c1ccc(OC)cc1` | 397.2 | 21 | MW=397.2 > 350.0 |
| `MQWLUYVHADZGGE-FEEDHTTJSA-N` | `N=C(O[C@H]1O[C@H](OCc2ccccc2)[C@@H](OCc2ccccc2)C(OCc2ccccc2)...` | 671.0 | 45 | MW=671.0 > 350.0 |
| `NBTYKHKRENOLBT-UHFFFAOYSA-N` | `COc1cc(I)c(OCC2C=Cc3ccccc3O2)cc1OC` | 424.2 | 23 | MW=424.2 > 350.0 |
| `NEUHVHWAWAEHDL-DAXYICLTSA-N` | `C[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](CO)O[C@H](Oc3cc(O)c4c(c...` | 580.5 | 41 | MW=580.5 > 350.0 |
| `NRJRMFXLQSNZGY-PWPASLGISA-N` | `O[C@H]1Cc2c(OCc3ccccc3)cc(OCc3ccccc3)cc2O[C@@H]1c1ccc(OCc2cc...` | 650.8 | 49 | MW=650.8 > 350.0 |
| `NXZCETLYSWWPGL-HZWQLJEPSA-N` | `O=C(O[C@H]1[C@H](O[C@H]2C(Br)O[C@H](COCc3ccccc3)[C@@H](OCc3c...` | 957.9 | 66 | MW=957.9 > 350.0 |
| `OKVSBVAPBXYLNJ-UHFFFAOYSA-N` | `COc1ccc(-c2cc(=O)c3c(OC)cc(OC)c(C4OC(CO)C(O)C(O)C4O)c3o2)cc1...` | 504.5 | 36 | MW=504.5 > 350.0 |
| `OOQYKGWPDLNECP-ZVZJHGLMSA-N` | `CC(=O)OC1O[C@H](C)C(OCc2ccccc2)[C@@H](OCc2ccccc2)[C@H]1OCc1c...` | 476.6 | 35 | MW=476.6 > 350.0 |
| `OZLPYHSAZFJYFY-HKTKGDTJSA-N` | `CC(=O)OC[C@@H]1O[C@@H](OCC(=O)c2ccc(OC(C)=O)cc2)[C@@H](OC(C)...` | 524.5 | 37 | MW=524.5 > 350.0 |
| `PHMOUQZDAZMQSN-UHFFFAOYSA-N` | `C[Si](C)(C)CCOc1ccccc1C#Cc1ccccc1Br` | 373.4 | 22 | MW=373.4 > 350.0 |
| `PVSQRRGIWVPVBW-UHFFFAOYSA-N` | `COc1ccc(-c2oc3cc(OCc4ccccc4)cc(O)c3c(=O)c2O)cc1O` | 406.4 | 30 | MW=406.4 > 350.0 |
| `QOWFBQWCVIRVON-UHFFFAOYSA-N` | `COc1ccc2c(=O)c(-c3ccc(OCc4ccccc4)cc3OCc3ccccc3)coc2c1` | 464.5 | 35 | MW=464.5 > 350.0 |
| `SGEWCQFRYRRZDC-UHFFFAOYSA-N` | `O=c1cc(-c2ccc(O)cc2)oc2c(C3OC(CO)C(O)C(O)C3O)c(O)cc(O)c12` | 432.4 | 31 | MW=432.4 > 350.0 |
| `SMRRYUGQTFYZGD-UHFFFAOYSA-K` | `CC(=O)O[Tl](OC(C)=O)OC(C)=O` | 381.5 | 13 | MW=381.5 > 350.0 |
| `SYUVAXDZVWPKSI-UHFFFAOYSA-N` | `CCCC[Sn](CCCC)(CCCC)c1ccccc1` | 367.2 | 19 | MW=367.2 > 350.0 |
| `UOJHCSLDVZMPAR-UHFFFAOYSA-N` | `COc1ccc2c(c1)OC(COc1ccccc1I)C=C2` | 394.2 | 21 | MW=394.2 > 350.0 |
| `UZLTVNHMAAJFGJ-OQCRIDGZSA-N` | `O=C(O[C@H]1O[C@H](Br)[C@H](OC(=O)c2ccccc2)[C@H](OC(=O)c2cccc...` | 645.5 | 43 | MW=645.5 > 350.0 |
| `WFDCFOGFWKJBMJ-GWEQUVLNSA-N` | `C[C@H]1O[C@H](Oc2cc(O)c3c(=O)c(O[C@H]4OC(CO)(CO)[C@H]4O)c(-c...` | 564.5 | 40 | MW=564.5 > 350.0 |
| `XBSYMIMBHPOPFD-XNTDXEJSSA-N` | `COc1ccc(/C=C/C(=O)c2ccc(OC)cc2OCc2ccccc2)cc1` | 374.4 | 28 | MW=374.4 > 350.0 |
| `XGEYXJDOVMEJNG-GKONUYDESA-N` | `O=C(/C=C/c1ccc(O)c(O)c1)c1ccc(O[C@H]2O[C@@H](CO)[C@H](O)[C@@...` | 450.4 | 32 | MW=450.4 > 350.0 |
