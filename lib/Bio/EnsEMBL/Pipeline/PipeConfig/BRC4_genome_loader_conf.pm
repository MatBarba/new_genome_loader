package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.4;
use Bio::EnsEMBL::ApiVersion qw(software_version);

use File::Spec::Functions qw(catdir);
use FindBin;

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    ##############################
    # Config to be set by the user
    # VB release id
    release => $self->o('release'),
    db_prefix => "",

    # Ensembl version (deduced from the environment?)
    ensembl_version => software_version(),

    # Meta configuration directory
    data_dir => $self->o('data_dir'),

    # Working directory
    tmp_dir => 'tmp',

    check_manifest => 1,
    prune_agp => 1,

    ##############################

    # Basic pipeline configuration
    pipeline_tag => '',
    pipeline_name => 'brc4_genome_loader_' .
      $self->o('release') . '_' . $self->o('ensembl_version') . $self->o('pipeline_tag'),
    email => $ENV{USER} . '@ebi.ac.uk',
    pipeline_dir => 'genome_loader_' .  $self->o('release') . '_' . $self->o('ensembl_version'),

    debug => 0,
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    debug          => $self->o('debug'),
    tmp_dir        => $self->o('tmp_dir'),
    check_manifest => $self->o('check_manifest'),
    pipeline_dir   => $self->o('pipeline_dir'),
    ensembl_root_dir => $self->o('ensembl_root_dir'),

    proddb_url   => $self->o('proddb_url'),
    taxonomy_url => $self->o('taxonomy_url'),
    dbsrv_url    => $self->o('dbsrv_url'),

    prune_agp    => $self->o('prune_agp'),
  };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('tmp_dir'),
    ];
}

# Ensures output parameters get propagated implicitly
sub hive_meta_table {
  my ($self) = @_;

  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return
  [
    # Starting point
    # Create the tmp pipeline directory
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'FillDBParams',
        'A->1' => 'Cleanup',
      },
    },

    {
      # Delete the temp working directory
      -logic_name => 'Cleanup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
    },

    {
      # fill db params from urls
      # output ..._host ..._port ...
      -logic_name => 'FillDBParams',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DbUrlToParams',
      -parameters => {
        db_urls => {
          dbsrv    => '#dbsrv_url#',
          proddb   => '#proddb_url#',
          taxonomy => '#taxonomy_url#',
        },
      },
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into => 'Manifest_factory',
    },

    {
      # Create a thread for each species = manifest file
      # Output:
      -logic_name        => 'Manifest_factory',
      -module         => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters        => {
        data_dir => $self->o('data_dir'),
        inputcmd => "find #data_dir# -type f -name manifest.json",
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -max_retry_count => 0,
      -flow_into => {
        2 => WHEN('#check_manifest#', {
          'Manifest_integrity' => { manifest => '#_0#' }
        }, ELSE({
          'Prepare_genome' => {manifest => '#_0#' }
        })),
      },
    },

    {
      # Check the integrity of the manifest before loading anything
      -logic_name => 'Manifest_integrity',
      -module     => 'Integrity',
      -language => 'python3',
      -analysis_capacity   => 5,
      -rc_name         => '8GB',
      -max_retry_count => 0,
      -flow_into => 'Prepare_genome',
    },

    {
      # Prepare all the metadata:
      # - species
      # - db_name
      # - manifest_metadata
      -logic_name => 'Prepare_genome',
      -module     => 'PrepareGenome',
      -language => 'python3',
      -parameters => {
        release => $self->o('release'),
        ensembl_version => $self->o('ensembl_version'),
        db_prefix => $self->o('db_prefix'),
      },
      -analysis_capacity   => 1,
      -max_retry_count => 0,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '2->A' => 'CleanUpAndCreateDB',
        'A->2' => 'LoadData',
      },
    },

    {
      # Init the Ensembl core
      -logic_name => 'CleanUpAndCreateDB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url'),
        sql     => [
          'DROP DATABASE IF EXISTS #db_name#;' ,
          'CREATE DATABASE #db_name#;' ,
        ],
      },
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into => [ 'LoadDBSchema' ],
    },

    {
      -logic_name => 'LoadDBSchema',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        input_file => $self->o('ensembl_root_dir') . '/ensembl/sql/table.sql',
      },
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into  => [ 'PopulateProductionTables' ],
    },

    {
      -logic_name    => "PopulateProductionTables",
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #dump_path#; ' .
          ' perl #base_dir#/ensembl-production/scripts/production_database/populate_production_db_tables.pl '
            . ' --host #dbsrv_host# --port #dbsrv_port# --user #dbsrv_user# --pass #dbsrv_pass# --database #db_name# '
            . ' --mhost #proddb_host# --mport #proddb_port# --muser #proddb_user# --mdatabase #proddb_dbname# '
            . ' --dumppath #dump_path# --dropbaks '
            . ' > #dump_path#/stdout 2> #dump_path#/stderr ',
        'base_dir'       => $self->o('ensembl_root_dir'),
        'dump_path' => $self->o('pipeline_dir') . '/#db_name#/create_core/fill_production_db_tables',
      },
    },

    {
      # Head analysis for the loading of data
      -logic_name => 'LoadData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'LoadSequenceData',
        'A->1' => 'LoadMetadata',
      },
    },

    {
      -logic_name => 'LoadSequenceData',
      -module     => 'LoadSequenceData',
      -language => 'python3',
      -parameters        => {
        work_dir => $self->o('pipeline_dir') . '/#db_name#/load_sequence',
        prune_agp => $self->o('prune_agp'),
      },
      -analysis_capacity   => 5,
      -rc_name         => '8GB',
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
      -flow_into  => [ 'PrepareAssemblyData' ],
    },

    {
      -logic_name => 'PrepareAssemblyData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'LoadAssemblyData',
      },
    },

    {
      -logic_name => 'LoadAssemblyData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'SetupAssemblyMetadata',
      },
    },

    {
      -logic_name => 'SetupAssemblyMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
    },

    {
      # Head analysis for the loading of the metadata
      -logic_name => 'LoadMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'FillMetadata',
        'A->1' => 'ConstructRepeatLib',
      },
    },

    {
      -logic_name => 'FillMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'FillTaxonomy',
      },
    },

    {
      -logic_name => 'FillTaxonomy',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'ConstructRepeatLib',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
    },
  ];
}

sub resource_classes {
    my $self = shift;
    return {
      'default'  	=> {'LSF' => '-q production-rh74 -M 4000   -R "rusage[mem=4000]"'},
      '8GB'       => {'LSF' => '-q production-rh74 -M 8000   -R "rusage[mem=8000]"'},
      '15GB'      => {'LSF' => '-q production-rh74 -M 15000  -R "rusage[mem=15000]"'},
      '32GB'  	 	=> {'LSF' => '-q production-rh74 -M 32000  -R "rusage[mem=32000]"'},
      '64GB'  	 	=> {'LSF' => '-q production-rh74 -M 64000  -R "rusage[mem=64000]"'},
      '128GB'  	 	=> {'LSF' => '-q production-rh74 -M 128000 -R "rusage[mem=128000]"'},
      '256GB'  	 	=> {'LSF' => '-q production-rh74 -M 256000 -R "rusage[mem=256000]"'},
	}
}

1;
