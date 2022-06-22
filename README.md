English md is in preparation.

分岐理論の詳細についてはaw02m/bifurcation_theoryを参照してください．(加筆中)

Aug/31/2021 : QT5-QCustomPlotを用いたppツールを追加しました．(pp_qt)  
Nov/4/2021 : 一時的にPD計算しか利用できません．変分方程式の計算を高速化しました(1000倍くらい)．  
Nov/4/2021 : PD,NSに対応，fixが統合されました．  
Nov/16/2021 : NSをBialternate積方式に変更，Gに対応，積分器にboostを選択可能．  
Jun/22/2022 : Qmakeを廃止し，Cmakeに完全に移行．現在LinuxとMacOSでの動作を確認しています．PythonバージョンのPPはどう頑張ってもとろこいのでリストラしました．**READMEはまだ更新しきれてませんので，記述がおかしいところがあります．**  

# autonomous_bif
自律系(微分方程式)の分岐解析ツールです．
本プログラムは2部に分かれています．

1. pp (Phase portrait; 相平面描画ツール)
2. bif (Bifurcation; 分岐集合の計算)

## pp_qt概要
pp_qtはppに比べ非常に高速な描画を提供します．

### 動作環境
* Eigen-3.4-rc1
* nlohmann-3.10.1
* qt6 (クロスプラットフォーム描画フレームワーク)

### pp_qt入力ファイル概要
ppの入力パラメタを引き継ぎますが，次の追加パラメタが使用できます．
* "direction" : Poincare断面を横切りを判定する際の向きを決定します．non-zero
* "period" : 指定した整数に従って周期τを求めます．
* "max_plot" : 解軌道を描画する最高打点数を指定します．
* "max_poincare_plot" : Poincare写像の最高打点数を指定します．
* "use_classic_rk" : trueの場合古典的Runge-Kutta法を使用します．falseの場合は刻み幅を自動計算するRunge-Kutta-Fehlberg法を用います．false推奨．
* "rk_div" : 古典的Runge-Kutta法を用いる場合の周期時刻分割数です．
* "rkf_first_h" : RKF法を用いる際の初期ステップを指定します．0.01~0.001を推奨．
* "rkf_h_max" : RKF法の最大ステップを指定します．ステップhはこれ以上大きくなりません．
* "rkf_h_min" : RKF法の最小ステップを指定します．ステップhはこれ以上小さくなりません．
* "rkf_tol" : RKF法の許容誤差(制度保証)を指定します．対象の系によってアドホックに変更する必要あり．
* "rkf_false_iter" : RKF法にて計算に失敗したと判定する繰り返し数を指定します．
* "poincare_eps" : 解軌道がPoincare断面を横切った際のポアンカレ写像の点はNewton法にて求められるため，その際の許容誤差を指定します.

### 使用法
```
cd pp
mkdir build
cd build
cmake ../cmake-tree
make
./main ../input/<your_input_file>.json
```
* ←→ : 変更するパラメタを選択します．
* ↑↓ : 選択したパラメタを変更します．
* Space : Poincare断面上の固定点，パラメタ，周期τを出力します．
* x : x軸に対応する変数を変更します．
* y : y軸に対応する変数を変更します．
* t : 解軌道の描画の削除/表示を切り替えます．
* p : Poincare写像の削除/表示を切り替えます．
* q : プログラムを終了します．

## bif概要
bifはppで取得した固定点情報をもとに，パラメタを変化させながらNewton法にて精度の良い固定点計算および分岐計算を行います．

### 動作環境
* Eigen-3.4-rc1 : 線形代数ライブラリです．比較的新しい関数を使用しているため，gitリポジトリの最新バージョンを利用してください．Arch Linuxなひとは`$ yay -S eigen-git`でインストールされます．
* nlohmann-3.10.1 : jsonライブラリ．
* cmake : Makefileの自動生成に用います．面倒な人はMakefileを自力で書いてください．

### fix入力ファイル概要
* "x0" : ppにてjsonファイルを出力した際に追加されます．この値が固定点計算の初期値として使用されます．
* "tau" : 周期解の周期時刻を入力します．
* "p_index" : ポアンカレ断面を配置する軸(変数)を指定します．
* "p_place" : ポアンカレ断面を配置する場所を指定します． q(x[p_index])-p_place = 0 に断面が配置されます．
* "sigma" : bif_nsにおける特性定数の実部を入力します．bifではPD:-1, G:1のどちらかを入力.
* "omega" : bif_nsにおける特性定数の虚部を入力します．bifでは使用しません．
* "var_param" : 変数パラメタを指定します．"params"のインデックスで指定してください．
* "inc_param" : 変化させるパラメタを指定します．"params"のインデックスで指定してください．
* "delta_inc" : パラメタ変分量です．経験的に0.1~0.001の範囲で設定すると良いかと思います．
* "inc_iter" : 指定した回数固定点を計算します．(指定した計算回数に到達もしくは発散しないとプログラムは固定点を計算し続けます．)
* "max_iter" : Newton法の最大ステップを指定します．Newton法は通常数回の繰り返しで収束するため，10~32の整数値を指定します．
* "eps" : Newton法の収束判定を与えます．誤差が指定数値以下になった場合に計算が終了します．
* "explode" : Newton法の発散判定を与えます．誤差が指定数値以上になった場合に計算が終了します．

### 使用法
```
cd bif
mkdir build
cd build
cmake ../cmake-tree
make
./main ../input/<your_input_file>.json
```
力学系の写像及びその微分は`sys_func.cpp`に記述してください．
固定点計算に成功すると固定点座標，パラメタ値，特性定数，特性定数のノルム・偏角が出力されます．特性定数のノルムが1に近いもの(分岐点)をピックアップしてbifプログラムに渡してください．