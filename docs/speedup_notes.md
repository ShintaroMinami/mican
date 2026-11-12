# 高速化メモ（ブランチ `perf/speedup`）

2026-11-12 開始。C コア（`src/`）の高速化に関する決定事項と計測結果を記録する。

## 方針

- **出力のビット一致を保つ**。浮動小数点の演算順序・タイブレーク・丸めは変えない。
  - 回帰チェックの内容:
    - 890 ケース（MALISAM `-w/-s/-r`、MALIDUP `-w/-s`、skill_test の大きい構造 `-w/-s/-r`）で、`-z` 標準出力を全文比較する。
    - `-o/-a/-m` の出力ファイル（124 個、チェーン指定ありを含む）を比較する。
    - arm64（FMA あり）と x86_64（Rosetta、FMA なし＝CI と同条件）の両方で、ベースラインと最適化版が一致することを確認する。
- アルゴリズムの結果を変える案は採用しない。採用しないもの:
  - `-ffast-math`、`expf/acosf`、`-march=native`（x86 で FMA が有効になる）
  - qsort の置き換え
  - 累積変数の double 化

## 実施した変更

| # | 対象 | 内容 |
|---|---|---|
| 1 | `main_align.c` | GH ハッシュテーブルを一括確保に変更（各セルの calloc 約 15 万回を、共有プールからの切り出しに置き換え） |
| 2 | `global_alignment.c`, `ca_alignment.c` | **疎行列化**。mat[Nite_loc] > 0 の contact セルだけをリストで保持する（詳細は下記） |
| 3 | `utils.c` ほか | 残基ごとの `float[4]` calloc を `calloc_coords()`（一括確保）に置き換え（TMscore, TMscore_mod, SPscore, rmsd, bstft） |
| 4 | `main_align.c` | 精密化ループの収束時打ち切り（詳細は下記） |
| 5 | `mican.h` | `distance_2` を `static inline` にした（全ペアから呼ばれるため） |
| 6 | `ca_alignment.c` | `angle()` の残基ごとのベクトルとノルムを事前計算する（演算は `angle()` と同一） |
| 7 | `superposition.c` | 作業用 RESDAT バッファを呼び出し間で再利用する |
| 8 | `global_alignment.c` | 前任セルが無い場合は番兵セル（index ncons）を指すようにし、calc_score の分岐を削減 |
| 9 | `read_pdb.c`, `fsscanf.c` | PDB の ATOM 行を固定カラムで直接読む。fsscanf と同じ意味論を保ち、数値は `strtof`/`strtol` で読む |
| 10 | `ca_alignment.c` | 行ごとの距離計算を SoA 配列上の分岐なしループにした（ベクトル化される） |

#2 疎行列化の詳細:
- 削除したもの: 密行列 4 枚と score/length の初期化、`check_close` の全走査、行・列単位の USED 塗り。
- 非 contact セルの値は常に 0 以下で、score は 0 になる。このため対角の前任セルが無い場合は「score 0・歩行停止」と等価になり、結果は厳密に一致する。

#4 収束時打ち切りの根拠:
- `ca_alignment` 後の `align_q` が反復前と同一なら、次の反復は同じ superposition を作る。
- その結果 `check_align` が失敗し、状態は現在のものに戻される。したがって次の反復は省略できる。

## 試して不採用にしたもの

- **三角不等式（主鎖の累積長）による距離計算のスキップ**: 遅くなったため不採用（大規模で 0.207 秒 → 0.311 秒）。sqrt を待つ逐次依存のコストが、元ループの並列実行（約 1ns/ペア）を上回った。
- **LTO（`-flto`）**: 3〜5% 速くなったが、ビルド環境への依存を避けるため不採用。代わりに #5 を実施した。

## 計測結果（Apple M4 Pro、CPU 時間、3 回の最小値）

| セット | ベースライン | 最適化後 | 倍率 |
|---|---|---|---|
| MALISAM+MALIDUP `-w`（371 件） | 8.95 s | 4.31 s | 2.1x |
| MALISAM `-s/-r`（260 件） | 4.91 s | 2.53 s | 1.9x |
| 大きい構造（18 件） | 7.10 s | 2.72 s | 2.6x |
| 1dfj × 5e29 `-s`（内部時間） | 0.806 s | 0.190 s | 4.2x |

`tests/test_malisam.py` と同じ判定（accuracy 0.6996、CSV 不一致 0）。

## 残りのホットスポット（最適化後）

- `ca_alignment` 約 46%。主に contact ペアでの acos/exp で、出力を変えずに削るのは困難。
- PDB 読込 約 10%。同じファイルを 3 回読んでいる。
- superposition/bstft 約 7%、global_alignment 約 6.5%。

## 今後の候補

- 候補単位（GH 候補・精密化）のスレッド並列化（pthreads）。結果は決定的に保てる。
  - 1 ジョブの壁時計時間を大きく短縮できる。
  - ただし pymican のバッチ実行ではプロセス並列と競合するため、オプトインにするのがよい。
- PDB を 1 回だけ読む構成への変更。
