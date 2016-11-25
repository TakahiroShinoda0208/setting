;;; ロードパスの追加
(setq load-path (append '("~/.emacs.d") load-path))

;;行数表示
(require 'linum)
(global-linum-mode)

;;; 対応する括弧を光らせる。
(show-paren-mode 1)

;;; diredから"r"でファイル名をインライン編集する
(require 'wdired)
(define-key dired-mode-map "r" 'wdired-change-to-wdired-mode)

;;カーソルがある関数名を表示させる
(which-function-mode 1)


;; c-modeやc++-modeなどcc-modeベースのモード共通の設定
(add-hook
 'c-mode-common-hook
 (lambda ()
   ;;関数の折りたたみを可能にする
   (hs-minor-mode 1)
   (define-key global-map (kbd "C-\\") 'hs-toggle-hiding)

   ;; BSDスタイルをベースにする
   (c-set-style "bsd")

   ;; スペースでインデントをする
   (setq indent-tabs-mode nil)

   ;; インデント幅を2にする
   (setq c-basic-offset 6)

   ;; 自動改行（auto-new-line）と
   ;; 連続する空白の一括削除（hungry-delete）を
   ;; 有効にする
   (c-toggle-auto-hungry-state 1)

   ;; CamelCaseの語でも単語単位に分解して編集する
   ;; GtkWindow         => Gtk Window
   ;; EmacsFrameClass   => Emacs Frame Class
   ;; NSGraphicsContext => NS Graphics Context
   (subword-mode 1)))


;; マウス・スクロールを滑らかにする（Mac Emacs 専用）
(setq mac-mouse-wheel-smooth-scroll t)

;; キーワードのカラー表示を有効化
;; 「t」の部分を「nil」にするとカラー表示をOffにできる
(global-font-lock-mode t)

;;特定のキーワードをハイライトにする(C-x w h :ハイライト、C-x w r :解除)
(global-hi-lock-mode 1)