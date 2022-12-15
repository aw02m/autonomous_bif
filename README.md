<!-- English md is in preparation. -->
Bifurcation analysis tool for autonomous dynamical systems.
Ask me if you need the English usage.

分岐理論の詳細についてはaw02m/bifurcation_theoryを参照してください．(加筆中)

# autonomous_bif
自律系(微分方程式)の分岐解析ツールです．
本プログラムは3部に分かれています．

1. cgen (c++ code generator; 系の方程式とその偏微分を出力)
2. pp (Phase portrait; 相平面描画ツール)
3. bif (Bifurcation; 分岐集合の計算)

## cgen概要
Python上でSympyを用いて自動的に系のC++コードを出力します．

### 動作環境
* Python3
* sympy

### cgenの設定
関数"func"内に系の方程式を記述します．状態変数は"x[i]", パラメタは"p[i]"です．
系の次元とパラメタの数を指定する"xdim, pdim"も指定してください．  
"python cgen.py"を実行すると自動的にppとbifに系とその偏微分が出力されます．

## pp概要
2次元相平面(phase-portrait)を描画します．

### 動作環境
* gcc/clang (LinuxとMacのみ確認)
* Eigen-3.4-rc1
* nlohmann-3.10.1
* qt6

### pp入力ファイル概要
* "json_out_path" : (プログラム終了時の)jsonデータの出力先
* "tick" : 1フレームに描画する時刻の範囲
* "max_plot" : 解軌道を描画する最高打点数
* "max_poincare_plot" : Poincare写像の最高打点数
* "axis" : x/y軸それぞれに表示する変数
* "period" : 指定した整数に従って周期τを計算
* "xrange", "yrange" : x/y軸の描画範囲
* "x0" : 初期値
* "params" : パラメタ
* "p_index", "p_place" : ポアンカレ断面をx[p_index] - p_place = 0に配置
* "direction" : Poincare断面を横切りを判定する際の向き 正負を実数で指定しますが，値は関係ないので-1とか1とかを適当に与えてください．
* "period" : 指定周期ごとに周期時刻を計算
* "use_classic_rk" : trueの場合古典的Runge-Kutta法．falseの場合はRunge-Kutta-Fehlberg法．false推奨．
* "rkf_first_h" : RKF法を用いる際の初期ステップ
* "rkf_h_max" : RKF法の最大ステップ ステップhはこれ以上大きくなりません．
* "rkf_h_min" : RKF法の最小ステップ ステップhはこれ以上小さくなりません．
* "rkf_tol" : RKF法の許容誤差(制度保証) 対象の系によってアドホックに変更する必要あり.
* "rkf_false_iter" : RKF法にて計算に失敗したと判定する繰り返し数
* "delta_inc" : プログラム実行中にインタラクティブにパラメタを変化させる量
* "max_iter" : Poincare断面上の固定点を計算する際に実行されるNewton法の最大繰り返し数
* "eps" : Newton法の許容誤差
* "explode" : Newton法の発散判定の上限値

### 使用法
```
cd pp
mkdir build
cd build
cmake ../cmake-tree
make
./main ../input/<your_input_file>.json
```
環境によっては build/main_autogen/include/ui_mainwindow.h のヘッダでコンパイルエラーが発生します．
その時は手動でそのファイルの18行目くらいを次のようにいじってください
```
#include "../../qcustomplot.h"
を
#include "qcustomplot.h"
に変更
```

* ←→ : 変更するパラメタを選択
* ↑↓ : 選択したパラメタを変更
* Space : Poincare断面上の固定点，パラメタ，周期τを出力
* x : x軸に対応する変数を変更
* y : y軸に対応する変数を変更
* s : 積分途中の最終状態を表示
* t : 解軌道の描画の削除/表示
* p : Poincare断面上固定点の削除/表示
* -+ : "delta_inc"を変化
* w : 固定点・パラメタ・周期の情報を更新してjsonを出力
* q : プログラムの終了

## bif概要
bifはppで取得した固定点/平衡点情報をもとに，パラメタを変化させながらNewton法にて精度の良い固定点計算および分岐計算を行います．

### 動作環境
* gcc/clang (LinuxとMacのみ確認)
* Eigen-3.4-rc1
* nlohmann-3.10.1
* cmake

### bif入力ファイル概要
* "mode" : 0-固定点 1-接線分岐(川上形式; dchi/dmu) 2-周期倍分岐 3-Neimark-Sacker分岐(Bialternate) 4-平衡点 5-平衡点の接線分岐 6-平衡点のHopf分岐
* "out_path" : データの出力先
* "json_out_path" : (プログラム終了時の)jsonデータの出力先
* "x0" : 固定点または平衡点の位置
* "tau" : 周期解の周期時刻(平衡点モードでは不使用)
* "params" : パラメタ
* "p_index", "p_place" : ポアンカレ断面を配置する場所を指定 q(x[p_index] - p_place)-p_place = 0
* "use_classic_rk" : 微分方程式と変分方程式のベクトル場の強さは極端に違う場合が多いので"true"を推奨
* "rkまわり" : pp参照
* "numerical_diff" : 第二変分を数値微分で計算
* "diff_strip" : "numerical_diff == true"の際の数値微分で用いる変化量
* "inc_param" : 変分パラメタを指定
* "var_param" : 変数パラメタを指定
* "delta_inc" : パラメタ変分量
* "inc_iter" : 指定した回数Newton法を計算
* "max_iter" : Newton法の最大反復数を指定
* "eps" : Newton法の収束判定
* "explode" : Newton法の発散判定

### 使用法
```
cd bif
mkdir build
cd build
cmake ../cmake-tree
make
./main ../input/<your_input_file>.json
```

### Change Log
Aug/31/2021 : QT5-QCustomPlotを用いたppツールを追加しました．(pp_qt)  
Nov/4/2021 : 一時的にPD計算しか利用できません．変分方程式の計算を高速化しました(1000倍くらい)．  
Nov/4/2021 : PD,NSに対応，fixが統合されました．  
Nov/16/2021 : NSをBialternate積方式に変更，Gに対応，積分器にboostを選択可能．  
Jun/22/2022 : Qmakeを廃止し，Cmakeに完全に移行．現在LinuxとMacOSでの動作を確認しています．PythonバージョンのPPはどう頑張ってもとろこいのでリストラしました．  
Jul/26/2022 : C++の関数記述を自動で行えるようにしました．cgen.pyにベクトル場を記述すればsys_func.cppは自動で生成されます．また，ppプログラムに解軌道の最後の状態を表示させる機能を追加しました．Sキーで表示されます．bifプログラムでは変分方程式の記述を変更したことにより大幅にパフォーマンスが向上しました．さらに本バージョンから二階変分の数値微分による計算に対応します．また平衡点計算が新たにbifに統合されました．(mode=4,5,6)  
Jul/30/2022 : また早くなりました．3次元系の周期解の分岐が10ms〜20msほどで求まります．PPプログラムのパラメタの変化量を-/+キーで調節できるようになりました．また，ポアンカレ断面上の固定点計算が条件によってスキップされるバグを修正しました．  
Jul/31/2022 : もっとはやくなりました．3次元系の周期解の分岐が1〜5msほどで求まります．変分方程式の記述を変更しました．どうでもいいですが妹が18歳になりました．かわいいでちゅね．  
Dec/12/2022 : 自律系の接線分岐の目的関数が不完全だった問題を解消しました．（自律系だと周期解計算ではμ=1が確定しているので，接線分岐条件が常に満たされてしまっていた．）具体的にはdchi/dmuを目的関数に用いる川上形式に変更．また，平衡点の接線分岐を実装するのを忘れていたので追加しました．
Dec/15/2022 : ppにjson出力機能を追加．また，ppとbifともに，jsonを出力した際にjsonオブジェクトが勝手にアルファベット順にソートされてしまう問題を修正しました．