# 少し凝った zshrc
# License : MIT
# http://mollifier.mit-license.org/

# aliasファイルを読み込む
if [ -f ~/.set/.alias_init ]; then
        . ~/.set/.alias_init
fi
# pathファイルを読み込む
if [ -f ~/.set/.path_init ]; then
        . ~/.set/.path_init
fi
# optionファイルを読み込む
if [ -f ~/.zsh/.option_init ]; then
        . ~/.zsh/.option_init
fi

#######################################
#基本コマンドのpath
export PATH=$PATH:$HOME/.local/bin:$HOME/bin

#言語の設定
LANG="en_US.UTF-8"



########################################
# 環境変数
#export LANG=ja_JP.UTF-8

# emacs 風キーバインドにする
bindkey -e

# 色を使用出来るようにする
autoload -Uz colors
colors

#colorの設定
if [ -f ~/.dircolors ]; then
    if type dircolors > /dev/null 2>&1; then
        eval $(dircolors ~/.dircolors)
    elif type gdircolors > /dev/null 2>&1; then
        eval $(gdircolors ~/.dircolors)
    fi
fi

# ヒストリの設定
HISTFILE=~/.zsh_history
HISTSIZE=1000000
SAVEHIST=1000000

# プロンプト
# 1行表示
# PROMPT="%~ %# "
#PROMPT="[%n%m](%*%) %~ %% "

# 2行表示
#PROMPT="%{${fg[white]}%}[%n@%m]%{${reset_color}%} %~
#%# "
PROMPT="%{${fg[white]}%}[%*%  %n@%m]%{${reset_color}%} %~
%# "




# 単語の区切り文字を指定する
autoload -Uz select-word-style
select-word-style default
# ここで指定した文字は単語区切りとみなされる
# / も区切りと扱うので、^W でディレクトリ１つ分を削除できる
zstyle ':zle:*' word-chars " /=;@:{},|"
zstyle ':zle:*' word-style unspecified

########################################
# 補完
# 補完機能を有効にする
autoload -Uz compinit
compinit -u

# 補完で小文字でも大文字にマッチさせる
zstyle ':completion:*' matcher-list 'm:{a-z}={A-Z}'

#補完時のファイルの色を表示する
if [ -n "$LS_COLORS" ]; then
    zstyle ':completion:*' list-colors ${(s.:.)LS_COLORS}
fi


# ../ の後は今いるディレクトリを補完しない
zstyle ':completion:*' ignore-parents parent pwd ..

# sudo の後ろでコマンド名を補完する
zstyle ':completion:*:sudo:*' command-path /usr/local/sbin /usr/local/bin \
                   /usr/sbin /usr/bin /sbin /bin /usr/X11R6/bin

# ps コマンドのプロセス名補完
zstyle ':completion:*:processes' command 'ps x -o pid,s,args'


########################################
# vcs_info
autoload -Uz vcs_info
autoload -Uz add-zsh-hook

zstyle ':vcs_info:*' formats '%F{green}(%s)-[%b]%f'
zstyle ':vcs_info:*' actionformats '%F{red}(%s)-[%b|%a]%f'

function _update_vcs_info_msg() {
    LANG=en_US.UTF-8 vcs_info
    RPROMPT="${vcs_info_msg_0_}"
}
add-zsh-hook precmd _update_vcs_info_msg

########################################
# OS 別の設定
case ${OSTYPE} in
    darwin*)
        #Mac用の設定
        #export CLICOLOR=1
	alias ls='gls --color'	
        ;;
    linux*)
        #Linux用の設定
           ;;
esac
