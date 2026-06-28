# Git & GitHub 工作流 SOP

> 本项目的版本控制规范与日常操作指南

---

## 1. 仓库信息

| 项目 | 值 |
|------|-----|
| 远程仓库 | https://github.com/juliusloon/retro_synthesis |
| 主分支 | `master` |
| 协议 | HTTPS / SSH |
| 用户 | cano (juliusloon) |

---

## 2. 每日工作流

### 2.1 开始工作

```bash
cd /home/ljx/retro_synthesis
git pull                    # 拉取远程最新代码
git status                  # 确认当前状态干净
```

### 2.2 工作中

```bash
# ... 编辑文件、运行脚本、做实验 ...

git status                  # 查看哪些文件变了
git diff                    # 查看具体修改内容
```

### 2.3 提交并推送

```bash
git add -A                  # 暂存所有修改
git commit -m "类型: 简要描述"  # 提交
git push                    # 推送到 GitHub
```

---

## 3. Commit 消息规范

### 3.1 格式

```
类型: 简要描述（不超过 50 字）

可选的详细说明（解释为什么做这个改动）
```

### 3.2 类型前缀

| 前缀 | 用途 | 示例 |
|------|------|------|
| `feat:` | 新功能、新模板、新脚本 | `feat: 添加异黄酮类反应模板 (45 条)` |
| `fix:` | 修复 bug | `fix: 修复 SMARTS 中的立体化学标记错误` |
| `docs:` | 文档更新 | `docs: 更新 SOP 中的配置文件路径` |
| `refactor:` | 重构代码（不改变功能） | `refactor: 合并 convert_*.py 为统一转换器` |
| `validate:` | 验证结果 | `validate: 橙皮苷模板匹配率从 78% 提升到 95%` |
| `data:` | 数据/模板文件更新 | `data: 更新 flavonoid_templates.hdf5 (508 → 540)` |
| `test:` | 测试配置/参数实验 | `test: 增加 MCTS 迭代到 1000，测试路线质量` |
| `chore:` | 杂项（环境、依赖等） | `chore: 更新 requirements.txt` |

### 3.3 好的 vs 不好的

```
# ✅ 好的
feat: 添加糖苷化反应模板，覆盖芦丁糖和新橙皮糖
fix: 修复 build_literature_templates.py 中 [R] 基团解析失败的问题
validate: 用 20 个代表性黄酮验证生物合成模板，全部通过

# ❌ 不好的
update
fix bug
提交
wip
```

---

## 4. 典型工作场景

### 4.1 添加新反应模板

```bash
# 1. 编辑学习笔记
vim files/26-06-28 新的黄酮反应.md

# 2. 运行提取脚本
python src/build_literature_templates.py

# 3. 验证模板
python src/validate_templates_on_hesperidin.py

# 4. 检查结果
git status
git diff --stat              # 只看文件变更统计

# 5. 提交
git add -A
git commit -m "feat: 从新笔记提取 30 条反应模板，橙皮苷验证通过"
git push
```

### 4.2 测试不同搜索参数

```bash
# 1. 备份当前配置
cp config/flavonoid.yml config/flavonoid.yml.bak

# 2. 修改参数
vim config/flavonoid.yml

# 3. 运行搜索
conda run -n retro aizynthcli \
  --smiles "COc1ccc(...)cc1O" \
  --config config/flavonoid.yml \
  --output outputs/26-06-28-test-1000iter.json

# 4. 查看结果，决定是否保留
python scripts/processing.py

# 5. 提交（如果结果好）
git add config/flavonoid.yml outputs/26-06-28-test-1000iter.json
git commit -m "test: MCTS 1000 迭代测试，top route 3 步"
git push

# 或者回退（如果结果不好）
mv config/flavonoid.yml.bak config/flavonoid.yml
```

### 4.3 开发新脚本

```bash
# 1. 创建新脚本
vim src/new_analysis.py

# 2. 测试
python src/new_analysis.py

# 3. 提交
git add src/new_analysis.py
git commit -m "feat: 添加新的分析脚本，用于 xxx"
git push
```

### 4.4 实验失败，回退修改

```bash
# 查看改了什么
git diff

# 放弃某个文件的修改（恢复到上次提交的状态）
git checkout -- config/flavonoid.yml

# 放弃所有修改
git checkout -- .

# 如果已经提交了但想回退
git log --oneline            # 找到要回退到的 commit
git revert HEAD              # 创建一个新 commit 来撤销上一次
# 或者
git reset --hard f55a21f     # 直接回到某个版本（危险！丢失后续修改）
```

---

## 5. 分支策略

### 5.1 什么时候用分支

| 场景 | 是否需要分支 |
|------|-------------|
| 小修改（修 bug、改配置） | 不需要，直接在 master 上提交 |
| 新增一批模板 | 不需要，直接提交 |
| 大规模重构（重写提取逻辑） | **建议用分支** |
| 实验性改动（不确定是否可行） | **建议用分支** |
| 多人协作同一功能 | **必须用分支** |

### 5.2 分支操作

```bash
# 创建并切换到新分支
git checkout -b feat/new-glycosylation-templates

# 在分支上工作、提交
git add -A
git commit -m "feat: 添加糖苷化模板"
git push -u origin feat/new-glycosylation-templates

# 完成后合并回 master
git checkout master
git pull                     # 先拉取 master 最新
git merge feat/new-glycosylation-templates
git push

# 删除分支（可选）
git branch -d feat/new-glycosylation-templates
git push origin --delete feat/new-glycosylation-templates
```

### 5.3 分支命名规范

```
feat/xxx    — 新功能
fix/xxx     — 修复
test/xxx    — 实验
refactor/xxx — 重构
```

---

## 6. `.gitignore` 说明

以下文件/目录 **不会** 被 Git 追踪:

| 类别 | 文件 |
|------|------|
| 大模型文件 | `data/*.onnx`, `data/*.hdf5`, `data/*.csv.gz` |
| 虚拟环境 | `venv/` |
| Python 缓存 | `__pycache__/`, `*.pyc` |
| 日志 | `*.log` |
| IDE | `.vscode/`, `.idea/` |
| OS | `.DS_Store` |

> ⚠️ `data/` 下的模型文件需要单独备份（约 1.7 GB）。建议保存到云盘或单独的存储位置。

---

## 7. 常用命令速查

### 查看状态

| 命令 | 作用 |
|------|------|
| `git status` | 当前状态（哪些文件变了） |
| `git diff` | 具体修改内容 |
| `git diff --stat` | 只看文件变更统计 |
| `git log --oneline` | 简洁历史 |
| `git log --oneline -10` | 最近 10 条 |
| `git log --oneline --graph` | 带分支图的历史 |
| `git remote -v` | 远程仓库地址 |

### 提交与推送

| 命令 | 作用 |
|------|------|
| `git add -A` | 暂存所有修改 |
| `git add <file>` | 暂存指定文件 |
| `git add -p` | 交互式选择要暂存的修改 |
| `git commit -m "msg"` | 提交 |
| `git commit --amend` | 修改上一次 commit |
| `git push` | 推送到远程 |

### 分支

| 命令 | 作用 |
|------|------|
| `git branch` | 列出本地分支 |
| `git branch -a` | 列出所有分支（含远程） |
| `git checkout -b <name>` | 创建并切换分支 |
| `git checkout master` | 切回 master |
| `git merge <branch>` | 合并分支 |

### 撤销与回退

| 命令 | 作用 |
|------|------|
| `git checkout -- <file>` | 放弃某个文件的修改 |
| `git stash` | 暂存当前修改 |
| `git stash pop` | 恢复暂存的修改 |
| `git revert HEAD` | 撤销上一次 commit（安全） |
| `git reset --hard <hash>` | 回退到指定版本（危险） |

---

## 8. 提交节奏建议

| 时机 | 说明 |
|------|------|
| 完成一个逻辑单元 | 提取完一批模板、修完一个 bug |
| 实验有结果后 | 搜索结果、验证报告值得记录 |
| 下班前 | 保存当天进度 |
| 大改动前 | 先 commit 当前状态，方便回退 |
| 跑长时间脚本前 | 确保代码已提交，避免丢失 |

**原则: 小步提交，频繁推送。**

---

## 9. 注意事项

1. **不要提交大文件**: `data/` 下的模型文件已被 `.gitignore` 排除，不要手动 `git add -f` 强制提交
2. **不要提交敏感信息**: API key、密码等不要写在代码里
3. **提交前先 `git diff`**: 确认修改内容是你预期的
4. **push 前先 pull**: 养成 `git pull` 再 `git push` 的习惯，避免冲突
5. **有意义的 commit 消息**: 未来的你会感谢现在写了清楚注释的自己
