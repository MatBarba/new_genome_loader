package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpGenomeJson;

use strict;
use warnings;

use JSON;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    metadata_name => 'genome',
  };
}

sub prepare_data {
  my ($self) = @_;

  my $sub_dir = $self->create_dir('json');
  $self->info(
          "Processing " . $self->production_name() . " genome into $sub_dir" );

  # Get meta table
  my $ma = Bio::EnsEMBL::Registry->get_adaptor($self->production_name, 'Core', 'MetaContainer');
  my $dba = $self->core_dba();

  my $meta = {
    species => {},
    assembly => {},
    genebuild => {}
  };

  my %meta_list = (
    species => [ qw(taxonomy_id production_name scientific_name strain display_name division alias) ],
    assembly => [ qw(accession date name version) ],
    genebuild => [ qw(version method start_date) ],
    provider => [ qw(name url) ],
  );
  for my $domain (keys %meta_list) {
    for my $key (@{$meta_list{$domain}}) {
      my @values = @{ $ma->list_value_by_key($domain . '.' . $key) };
      next if scalar(@values) < 1;
      @values = map { ($_ =~ /^\d+$/)? int($_) : $_ } @values;
      $meta->{$domain}->{$key} = (scalar(@values) > 1)? \@values : $values[0];
    }
  }

  # Check the assembly version
  $meta = $self->check_assembly_version($meta);

  $dba->dbc()->disconnect_if_idle();

  return $meta;
}

sub check_assembly_version {
  my ($self, $meta) = @_;

  my $assembly = $meta->{assembly};

  # Is there a version?
  my $version = $assembly->{version} // $assembly->{name};
  if ($version) {
    # Version is an integer
    if (not $version =~ /\D/) {
      $meta->{assembly}->{version} = int($meta->{assembly}->{version});
      return $meta;
    }

    # Version is not an integer, but ends in one
    if ($version =~ /^[A-z]+([0-9]+)$/) {
      $meta->{assembly}->{version} = int($1);
      return $meta;
    
    # There is a version but I can't get a number out of it
    } else {
      die("Can't extract version number from assembly version: $version");
    }
  } else {
    use Data::Dumper;
    die("No assembly version found in meta: " . Dumper($meta));
  }
  return;
}

1;
