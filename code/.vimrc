set noswapfile
set autoindent
set cindent
set shiftwidth=4
set relativenumber
set number
set numberwidth=5
"Things that make python programing easier
syntax on
filetype indent plugin on
set modeline
set background=dark
set ruler	

"ABBREVIATIONS
iabbrev shp #!/usr/bin/env python3
iabbrev shb #!/bin/bash
iabbrev shr #!/usr/bin/env R
iabbrev ini #!/usr/bin/env python3<cr><cr>__appname__ = '[App_name_here]'<cr>__author__ = 'Pablo Lechon (plechon@ucm.es)'<cr>__version__ = '0.0.1'<cr><cr>## IMPORTS ##<cr><cr>import sys<cr><cr>## CONSTANTS ##<cr><cr><cr>## FUNCTIONS ##<cr><cr>def main(argv):<cr>'''Main function'''<cr><cr>return 0<cr><cr>## CODE ##<cr><cr>if (__name__ == '__main__'):<cr>status = main(sys.argv)<cr>sys.exit(status)<cr> 
iabbrev sqares squares
iabbrev sqare square
iabbrev inic #include <stdio.h><cr>#include <stdlib.h><cr><cr>int main (void)<cr>{<cr><cr>return 0;<cr>}
iabbrev null NULL

autocmd FileType python :iabbrev <buffer> iff if :<left>
autocmd FileType r :iabbrev <buffer> iff if(){}<left><left><left>
"Leaders
let mapleader = ","
let localleader = "\\"

"MAPPINGS"
"Map the esc key to no operation to change habits
:inoremap <esc> <nop>
:inoremap kj <esc>
"Movement
"End and begining of current line
nnoremap <leader>L $a
nnoremap <leader>H 0i

"To move the current line up and down
noremap <leader>_ ddp
noremap <leader>- ddkkp

"Capitalize current word in normal and insert mode
nnoremap <leader><c-u> viwU
inoremap <leader><c-u> <esc>viwUea

"Inmediately entering and sourcing into vimrcfile
nnoremap <leader>ev :vsplit $MYVIMRC<cr>G
nnoremap <leader>sv :source $MYVIMRC<cr>

"Surround the cursor word with things
nnoremap <leader>" viw<esc>a"<esc>hbi"<esc>lel
nnoremap <leader>' viw<esc>a'<esc>hbi'<esc>lel
nnoremap <leader>( viw<esc>a)<esc>hbi(<esc>lel
nnoremap <leader>{ viw<esc>a}<esc>hbi{<esc>lel
nnoremap <leader>[ viw<esc>a]<esc>hbi[<esc>lel
nnoremap <leader>$ viw<esc>a$<esc>hbi$<esc>lel

"Surround the visually selected word/region with things
vnoremap <leader>" <esc>a"<esc>`<i"<esc>`>
vnoremap <leader>) <esc>a)<esc>`<i(<esc>`>
vnoremap <leader>[ <esc>a]<esc>`<i[<esc>`>
vnoremap <leader>{ <esc>a}<esc>`<i{<esc>`>
vnoremap <leader>$ <esc>a$<esc>`<i$<esc>`>
vnoremap <leader>' <esc>a'<esc>`<i'<esc>`>

"comment current line for different filetypes
autocmd FileType tex,plaintex nnoremap <leader>c I%<esc>
autocmd FileType python,r,sh nnoremap <leader>c I#<esc>
"comment multiple lines
vnoremap <leader>c i#<esc>

"insert blank line above/below and keep writing
nnoremap <leader>o o<esc>ja
nnoremap <leader>o o<esc>ka

"run .py, .r, .sh, .c, and .tex files
autocmd FileType r nnoremap <leader>r :w<cr>:!rscript '%'<cr>
autocmd FileType python nnoremap <leader>r :w<cr>:!python3 '%'<cr>
autocmd FileType sh nnoremap <leader>r :w<cr>:!bash '%'<cr>
"%:r prints the name of the file witout extension
autocmd FileType tex nnoremap <leader>r :w<cr>:w ! bash ~/Desktop/Imperial/CMEECourseWork/Week1/Code/CompileLaTex.sh '%:r'<cr>
"run a python program with arguments passed through the terminal
nnoremap <leader>a :w<cr>:!python3 '%' 

"change buffer quicker
nnoremap <leader>b :w<cr>:ls<cr>:b<space>
"change to previous buffer
nnoremap <leader>bb :w<cr>:ls<cr>:b<space>#<cr><cr>

"move between windows
nnoremap <leader>j <c-w>j
nnoremap <leader>k <c-w>k
nnoremap <leader>h <c-w>h
nnoremap <leader>h <c-w>h
nnoremap <leader>l <c-w>l
nnoremap <leader>l <c-w>l

"switch between normal mode and the terminal
tnoremap <esc> <c-w>n
nnoremap <leader>te <c-\>i 

"find and replace quickly
nnoremap <leader>fr :%s///gc<left><left><left><left>

"word count in latex
nnoremap <leader>wc :w !detex \| wc -w<cr>

"Insert name of current file
nnoremap <leader>cf "%p

"Insert date at the end of page with format
nnoremap <leader>d Go<CR>#### <C-R>=strftime("%b %d,%Y")<CR><Esc>o<CR>

"Git commands
nnoremap <leader>st :w ! git status<cr>
nnoremap <leader>am :! git add -A<cr> :! git commit -m ''<left>
nnoremap <leader>pu :w<cr>:! git push -u origin main<cr><cr>
nnoremap <leader>pd :w<cr>:! git add -A<cr>:! git commit -m '.'<cr>:! git push origin master<cr><cr>

"Save quickly
nnoremap <leader>w :w<cr>
inoremap <leader>w <esc>:w<cr>


"Insert a trace quickly"
nnoremap <S-d> oimport ipdb; ipdb.set_trace(context = 20)<Esc>
"Bring bottom line next to current line"
nnoremap <leader>u j0whvk$d
"Draw a gray margin line"
augroup python
    au!
    autocmd FileType python set colorcolumn=80 
    highlight ColorColumn ctermbg=lightgrey guibg=lightgrey
augroup END

augroup c
    au!
    autocmd FileType c set colorcolumn=80 
    highlight ColorColumn ctermbg=lightgrey guibg=lightgrey
augroup END
