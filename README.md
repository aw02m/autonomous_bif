English md is in preparation.
References are at of the end.

分岐理論の詳細についてはaw02m/bifurcation_theoryを参照してください．(加筆中)

# autonomous_bif
自律系(微分方程式)の分岐解析ツールです．
本プログラムは3部に分かれています．

1. pp (Phase portrait; 相平面描画ツール)
2. fix (Fixed point; 固定点計算及び特性定数の表示)
3. bif (Bifurcation; 分岐集合の計算)

## pp概要
ppは相平面をリアルタイムに描画します．(Original:Tetsushi Ueta)
固定点計算のための近似値取得に使用してください．

### 動作環境
* numpy
* scipy
* matplotlib

### pp入力ファイル要素
* "tick" : 一回の描画で進める積分時間
* "axis" : 横軸，縦軸に表示する変数
* "xrange", "yrange" : x,y描画区間．
* "x0" : プログラム実行時の初期値．
* "params" : 力学系のパラメタ．
* "dparams" : プログラム内でパラメタを変更する際の刻み幅．
* "p_index" : ポアンカレ断面を配置する軸(変数)を指定します．
* "p_place" : ポアンカレ断面を配置する場所を指定します．
** q(x[p_index])-p_place = 0 に断面が配置されます．

### 使用法
はじめに，pp.py内のシステム方程式を設定します．
その後，`$ python pp.py [json input file]`で実行．プログラム実行中は次のキーを押下することで各種パラメタ変更，ファイル出力等が行えます．
* `w` : `__ppout__.json`ファイルに現在の情報を保存します．ファイルがすでに存在している場合は上書きされます．
* `p` : 変化させるパラメタを指定します．現在の指定インデックスはターミナルに表示されます．
* `[矢印キー]` : 指定しているパラメタを"dparams"に従って変化させます．
* `[space]` : 描画をリセットします．
* その他の機能についてはProf. Uetaのオリジナルリポジトリを参照ください．

## fix概要
fixはppで取得した固定点情報をもとに，パラメタを変化させながらNewton法にて精度の良い固定点計算を行います．

### 動作環境
* Eigen3 : 線形代数ライブラリです．比較的新しい関数を使用しているため，gitリポジトリの最新バージョンを利用してください．Arch Linuxなひとは`$ yay -S eigen-git`でインストールされます．
* nlohmann : jsonライブラリ．
* cmake : Makefileの自動生成に用います．面倒な人はMakefileを自力で書いてください．

### fix入力ファイル概要
* "fixed" : ppにてjsonファイルを出力した際に追加されます．この値が固定点計算の初期値として使用されます．
* "tau" : 周期解の周期時刻を指定します．
* "p_index" : ポアンカレ断面を配置する軸(変数)を指定します．
* "p_place" : ポアンカレ断面を配置する場所を指定します．
** q(x[p_index])-p_place = 0 に断面が配置されます．
* "use_classic_rk" : trueの場合古典的Runge-Kutta法をもちいます．falseの場合は刻み幅を自動計算するRunge-Kutta-Fehlberg法を用います．false推奨．
* "rk_div" : 古典的Runge-Kutta法を用いる場合の周期時刻分割数です．
* "rkf_first_h" : RKF法を用いる際の初期ステップを指定します．0.01~0.001を推奨．
* "rkf_h_max" : RKF法の最大ステップを指定します．ステップhはこれ以上大きくなりません．
* "rkf_h_min" : RKF法の最小ステップを指定します．ステップhはこれ以上小さくなりません．
* "rkf_tol" : RKF法の許容誤差(制度保証)を指定します．対象の系によってアドホックに変更する必要あり．
* "rkf_false_iter" : RKF法にて計算に失敗したと判定する繰り返し数を指定します．
* "inc_param" : 変化させるパラメタを指定します．"params"のインデックスで指定してください．
* "delta_inc" : パラメタ変分量です．経験的に0.1~0.001の範囲で設定すると良いかと思います．
* "inc_iter" : 指定した回数固定点を計算します．(指定した計算回数に到達もしくは発散しないとプログラムは固定点を計算し続けます．)
* "max_iter" : Newton法の最大ステップを指定します．Newton法は通常数回の繰り返しで収束するため，10~32の整数値を指定します．
* "eps" : Newton法の収束判定を与えます．誤差が指定数値以下になった場合に計算が終了します．
* "explode" : Newton法の発散判定を与えます．誤差が指定数値以上になった場合に計算が終了します．

### 使用法
`./main [input json file]`にて実行．コンパイルは`fix`ディレクトリに新たに`build`ディレクトリを作成し，その中で`cmake`するとMakefileが自動で作成されます．  
力学系の写像及びその微分は`ds_derivatives.cpp`に記述してください．
計算に成功すると固定点座標，パラメタ値，特性定数，特性定数のノルム・偏角が出力されます．特性定数のノルムが1に近いもの(分岐点)をピックアップしてbifプログラムに渡してください．

## bif概要
完成しました．READMEは準備中...

# References
[1] Yuri, Kuznetsov, "Elements of Applied Bifurcation Theory", Springer, AMS 112.  
[2] T. Kousaka, T. Ueta, and H. Kawakami, “Bifurcation of switched nonlinear dynamical systems,” IEEE Transactions, CASII, vol. 46, no. 7, pp. 878–885, 1999.  
[3] T. UETA, M. TSUEIKE, H. KAWAKAMI, T. YOSHINAGA, and Y. KATSUTA, “A computation of bifurcation parameter values for limit cycles,” IEICE Transactions on Fundamentals of Electronics, Communications and Computer Sciences, vol. 80, no. 9, pp. 1725–1728, 1997, issn: 09168508.  
[4] T. Ueta, T. Yoshinaga, H. Kawakami, and G. Chen, “A method to calculate neimark-sacker bifurcation in autonomous systems,” IEICE-Transactions Fundamentals, EA, vol. 83, no. 10, pp. 1141–1147, 2000.  
[5] M. Bernardo, C. Budd, A. R. Champneys, and P. Kowalczyk, Piecewise-smooth dynamical systems: theory and applications. 2008, vol. 163.  
[6] S. Amoh, "A Bifurcation Analysis Method for High-dimensional Hybrid Autonomous Systems", Master thesis, Tokushima University, 2021.  
[7] H. Kawakami, "非線形現象 - 定常的接近法", 2005.
