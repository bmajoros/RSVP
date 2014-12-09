package EssexToken;
use strict;

######################################################################
#
# EssexToken.pm bmajoros@duke.edu 10/28/2008
#
# Represents a token in a Essex file --- either a parenthesis or a
# literal (a literal is anything other than a parenthesis: a tag/ident-
# ifier, number, or string).  Note that the scanner can't differentiate
# between a tag and a datum, so the parser must do this.
#
# Attributes:
#   tokenType : one of "(", ")", or "L" (literal)
#   lexeme : string
# Methods:
#   $token=new EssexToken($type,$lexeme);
#   $type=$token->getType();
#   $string=$token->getLexeme();
#   $bool=$token->isOpenParen();
#   $bool=$token->isCloseParen();
#   $bool=$token->isLiteral();
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $token=new EssexToken($type,$lexeme);
sub new
{
  my ($class,$type,$lexeme)=@_;

  my $self=
    {
     type=>$type
    };
  bless $self,$class;
  if(defined($lexeme) && length($lexeme)>0) {
    $self->{lexeme}=$lexeme
  }

  return $self;
}
#---------------------------------------------------------------------
#   $type=$token->getType();
sub getType
  {
    my ($self)=@_;
    return $self->{type};
  }
#---------------------------------------------------------------------
#   $string=$token->getLexeme();
sub getLexeme
  {
    my ($self)=@_;
    if(defined($self->{lexeme})) {
      return $self->{lexeme}
    }
  }
#---------------------------------------------------------------------
#   $bool=$token->isOpenParen();
sub isOpenParen
  {
    my ($self)=@_;
    return $self->{type} eq "(";
  }
#---------------------------------------------------------------------
#   $bool=$token->isCloseParen();
sub isCloseParen
  {
    my ($self)=@_;
    return $self->{type} eq ")";
  }
#---------------------------------------------------------------------
#   $bool=$token->isLiteral();
sub isLiteral
  {
    my ($self)=@_;
    return $self->{type} eq "L";
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

