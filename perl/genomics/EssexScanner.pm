package EssexScanner;
use strict;
use EssexToken;

######################################################################
#
# EssexScanner.pm bmajoros@duke.edu 10/28/2008
#
# A token scanner for the EssexParser.
# 
#
# Attributes:
#   file : FileHandle
#   ungot : char
#   nextToken : EssexToken
# Methods:
#   $scanner=new EssexScanner($filehandle);
#   $scanner->close();
#   $bool=$scanner->eof();
#   $token=$scanner->nextToken();
#   $token=$scanner->match($tokenType);
#   $token=$scanner->peek();
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $scanner=new EssexScanner($filehandle);
sub new
{
  my ($class,$file)=@_;

  my $self=
    {
     file=>$file
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $scanner->close();
sub close
  {
    my ($self)=@_;
    close($self->{file});
  }
#---------------------------------------------------------------------
#   $bool=$scanner->eof();
sub eof
  {
    my ($self)=@_;
    return eof($self->{file});
  }
#---------------------------------------------------------------------
#   $token=$scanner->nextToken();
sub nextToken
  {
    my ($self)=@_;
    my $token=$self->peek();
    undef($self->{nextToken});
    return $token;
  }
#---------------------------------------------------------------------
#   $token=$scanner->match($tokenType);
sub match
  {
    my ($self,$tokenType)=@_;
    my $token=$self->nextToken();
    if($token->getType() ne $tokenType) {
      my $lexeme=$token->getLexeme();
      die "Syntax error near \"$lexeme\"\n";
    }
  }
#---------------------------------------------------------------------
#   $token=$scanner->peek();
sub peek
  {
    my ($self)=@_;
    if(!defined($self->{nextToken})) {
      my $file=$self->{file};
      $self->skipWhitespace();
      my $c=$self->getChar();
      if(!defined($c)) { return undef }
      my ($tokenType,$lexeme);
      if($c=~/[\(\)]/) { $tokenType=$c }
      else {
	$tokenType="L";
	$lexeme=$c;
	while(!eof($file)) {
	  $c=$self->getChar();
	  if($c=~/[\s\(\)]/) {last}
	  $lexeme.=$c;
	}
	$self->unGetChar($c);
      }
      $lexeme=~s/&lparen;/\(/g;
      $lexeme=~s/&rparen;/\)/g;
      $lexeme=~s/&tab;/\t/g;
      $lexeme=~s/&space;/ /g;
      $self->{nextToken}=new EssexToken($tokenType,$lexeme);
    }
    return $self->{nextToken};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
sub getChar
  {
    my ($self)=@_;
    my $c;
    if(defined($self->{ungot})) {
      $c=$self->{ungot};
      undef($self->{ungot});
    }
    else {
      $c=getc($self->{file});
    }
    return $c;
  }
#---------------------------------------------------------------------
sub unGetChar
  {
    my ($self,$c)=@_;
    $self->{ungot}=$c;
  }
#---------------------------------------------------------------------
sub skipWhitespace
  {
    my ($self)=@_;
    my $file=$self->{file};
  WS:
    while(!eof($file)) {
      my $c=$self->getChar();
      last unless defined($c);
      if($c eq "#") {
	while(!eof($file)) {
	  $c=$self->getChar();
	  if(!defined($c) || $c=~/\n/) {last}
	}
	next WS;
      }
      if($c!~/\s/) {
	$self->unGetChar($c);
	return;
      }
    }
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------


1;

