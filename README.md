# Tools for NINJA MC study

NINJA 実験物理ランA (J-PARC E71a) における MC スタディのためのプログラム
以下のソフトウェアが使えるようになっていることが前提

 - GEANT4 v10.5.0+
 - ROOT 6.20+
 - BOOST 1.53+
 - [WAGASCI/Baby MIND Monte Carlo](https://git.t2k.org/wagasci_babymind/wagasci-babymind-monte-carlo) v0.1.12+
 - [WAGASCI/Baby MIND event display](https://git.t2k.org/wagasci_babymind/wagasci-babymind-event-display) v0.1.0+
 - [WAGASCI/Baby MIND track reconstruction](https://git.t2k.org/wagasci_babymind/wagasci-babymind-track-reconstruction) v0.0.2+
詳しい設定方法は WAGASCI/Baby MIND Monte Carlo のREADME.md を参照のこと

 - [NINJA track reconstruction tool](https://github.com/t-odagawa/ninjarecon)
詳しい設定方法は [NINJA wiki](https://www-he.scphys.kyoto-u.ac.jp/research/Neutrino/ninja/dokuwiki/doku.php?id=start) などを参照のこと

## 簡単な説明

### TrueDistribution
検出器アクセプタンスを考慮しない NEUT true base の分布を書くためのソフトウェア

### MatchAcceptance
検出器アクセプタンスを study するためのソフトウェア，シフター・トラッカーペーパーの図の作成に使用した

### Efficiency
検出効率を study するためのソフトウェア，データを用いて求めた efficiency を評価する？
[NinjaTrackerPerformance](https://github.com/t-odagawa/ninjatrackerperformance) も参照のこと

### EventSelection
Minimum Distance を始めとした Event selection を study するためのソフトウェア

### Kink
特に kink event を study するためのソフトウェア

### Background
Beam related な background を study するためのソフトウェア

### ReconDistribution
検出器アクセプタンスや検出効率を考慮した reconstructed base の分布を書くためのソフトウェア

