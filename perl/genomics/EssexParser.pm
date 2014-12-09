package EssexParser;
use strict;
use FileHandle;
use EssexNode;
use EssexScanner;
use EssexToken;

######################################################################
#
# EssexParser.pm bmajoros@duke.edu 10/27/2008
#
# A parser for the Essex language, which is a simple alternative to XML
# and looks like LISP programs ("S-Expressions" -- hence the name "Essex").
# It generates tree data structures built out of EssexNode objects.
#
# Attributes:
#   file : FileHandle
#   isOpen : boolean
#   scanner : EssexScanner
# Methods:
#   $parser=new EssexParser($filename); # filename is optional
#   $parser->open($filename); # unnecessary if you gave filename to ctor
#   $parser->close();
#   $tree=$parser->nextElem();   # returns root of the tree
#   $forest=$parser->parseAll(); # returns an array of trees
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class,$filename)=@_;

  my $self=
    {
     isOpen=>0
    };
  bless $self,$class;
  if($filename) { $self->open($filename) }

  return $self;
}
#---------------------------------------------------------------------
#   $parser->open($filename);
sub open
  {
    my ($self,$filename)=@_;
    if($self->{isOpen}) {$self->close()}
    my $file=$self->{file}=new FileHandle($filename);
    $self->{isOpen}=1;
    $self->{scanner}=new EssexScanner($file);
  }
#---------------------------------------------------------------------
#   $parser->close();
sub close
  {
    my ($self)=@_;
    if($self->{isOpen}) {
      close($self->{file});
      $self->{isOpen}=0;
      undef($self->{scanner});
    }
  }
#---------------------------------------------------------------------
#   $forest=$parser->parseAll(); # returns an array of trees
sub parseAll
  {
    my ($self)=@_;
    my $file=$self->{file};
    my $forest=[];
    while(!eof($file)) {
      my $tree=$self->nextElem();
      push @$forest,$tree;
    }
    return $forest;
  }
#---------------------------------------------------------------------
#   $tree=$parser->nextElem();   # returns root of the tree
sub nextElem
  {
    my ($self)=@_;
    die "file is not open\n" unless $self->{isOpen};
    my $scanner=$self->{scanner};
    my $token=$scanner->nextToken();
    if(!defined($token)) { return undef }
    if($token->isOpenParen()) { return $self->parseTuple() }
    elsif($token->isLiteral()) { return $token->getLexeme() }
    else {
      my $lexeme=$token->getLexeme();
      die "Syntax error near \n$token\"\n";
    }
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

sub parseTuple
  {
    # PRECONDITION: a "(" has already been matched

    my ($self)=@_;
    my $file=$self->{file};
    my @elements;
    my $scanner=$self->{scanner};
    while(!eof($file)) {
      my $token=$scanner->nextToken();
      #print $token->getType(); print "\n";
      last unless defined($token);
      if($token->isOpenParen()) {
	push @elements,$self->parseTuple();
      }
      elsif($token->isCloseParen()) { last }
      else {
	push @elements,$token->getLexeme();
      }
    }
    return new EssexNode(@elements);
  }


1;

