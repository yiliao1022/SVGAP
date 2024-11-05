package Parallel::ForkManager::Child;
our $AUTHORITY = 'cpan:DLUX';
# ABSTRACT: role adopted by forked Parallel::ForkManager processes 
$Parallel::ForkManager::Child::VERSION = '2.02';
    

use strict;
use warnings;

use Carp;

use Moo::Role;

sub is_child  { 1 }
sub is_parent { 0 }

sub start {
    croak "Cannot start another process while you are in the child process";
}

sub finish {
  my ($s, $x, $r)=@_;

  $s->store($r);

    CORE::exit($x || 0);
}

1;

__END__

=pod

=encoding UTF-8

=head1 NAME

Parallel::ForkManager::Child - role adopted by forked Parallel::ForkManager processes 

=head1 VERSION

version 2.02

=head1 SYNOPSIS 

   use 5.10.0;

   use Parallel::ForkManager;

   my $fm = Parallel::ForkManager->new;

   say "parent does not consume the child role: ", $fm->does('Parallel::ForkManager::Child');

   $fm->start_child(sub{ 
        sleep $_;
        say "but the child does: ", $fm->does('Parallel::ForkManager::Child');
        say "child $_ says hi!"

   }) for 1..3;

=head1 DESCRIPTION 

When the parent L<Parallel::ForkManager> object forks a child process,
its forked incarnation consumes this role. The role doesn't do much: it 
changes the returning values of C<is_child> and C<is_parent> in the way you'd expect,
change C<start> so that it'd die if called from within the child, and change 
the implementation of C<finish> to potentially send data back to the parent process.

=head1 AUTHORS

=over 4

=item *

dLux (Szabó, Balázs) <dlux@dlux.hu>

=item *

Yanick Champoux <yanick@cpan.org>

=item *

Gabor Szabo <gabor@szabgab.com>

=back

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2018, 2016, 2015 by Balázs Szabó.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
