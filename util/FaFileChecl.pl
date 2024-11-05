#!/usr/bin/perl

use strict;
use warnings;

# 从命令行参数获取目录路径
my $dir = $ARGV[0] or die "Usage: $0 <directory path>\n";

# 打开目录句柄
opendir(my $dh, $dir) or die "Cannot open directory: $!";

# 读取目录下的所有文件（排除以点开头的隐藏文件）
my @files = grep { !/^\./ && -f "$dir/$_" } readdir($dh);

# 关闭目录句柄
closedir($dh);

foreach my $file (@files) {
    # 检查文件名是否不以数字开头
    if ($file =~ /^\d/) {
        print "File starts with a number, which is not allowed: $file\n";
        next;
    }

    # 打开并读取文件，检查是否符合FASTA格式且序列ID匹配文件名
    open(my $fh, '<', "$dir/$file") or die "Cannot open file: $!";
    my $is_fasta = 0;
    my $id_matches_filename = 0;
    my $fasta_id;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)/) {
            $is_fasta = 1;
            $fasta_id = $1;
            last; # 只检查第一个FASTA头行
        }
    }
    close($fh);

    unless ($is_fasta) {
        print "File does not appear to be in FASTA format: $file\n";
        next;
    }

    # 提取文件名部分，不包括扩展名
    # my ($filename_without_ext) = $file =~ /^(.*)\.[^.]+$/;
    # unless ($fasta_id eq $filename_without_ext){
	# my ($filename_without_ext) = $file;
    unless ($fasta_id =~/$file/) {
        print "File name does not match FASTA ID: $file\n";
    }
}

print "Finished checking files.\n";

