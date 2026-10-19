***************************************************************************
*
* Title:                 ChEBI flat files
*
* ChEBI Release:         256
*
* Date of last update:   2026-10-06
*
***************************************************************************

The data content in ChEBI is licensed under a highly permissive Creative
Commons license - specifically the "CC Attribution-ShareAlike 4.0 International License.",
see LICENSE file. The required attribution should contain the url
of the ChEBI resource, and also the release version, e.g.:

ChEBI data is from http://www.ebi.ac.uk/chebi - the version of ChEBI is
256.

This should be visible on the entry portal for a web resource in which
ChEBI is integrated, or contained with the documentation for any further
distribution.

In this folder, you will find the same tables from the SQL dumps, but in TSV (tab-separated values) format,
with a slight difference: the TSV files contain only the primary IDs. As you may know, a compound in ChEBI can have
secondary IDs that reference the primary ID. You can find more information about this, as well as detailed documentation
on the columns here: https://chembl.gitdocs.ebi.ac.uk/chebi/chebi-2.0/chebi-data-management/#!/overview (Work In Progress)

Contents of this directory:
---------------------------

README                               This file

LICENSE                              License data is released under:
                                     Creative Commons Attribution 4.0 International License.

compounds.tsv.gz                     A compressed file of the compounds table separated by tabs (\t)

chemical_data.tsv.gz                 A compressed file of the chemical_data table separated by tabs (\t)

compound_origins.tsv.gz              A compressed file of the compound_origins table separated by tabs (\t)

structures.tsv.gz                    A compressed file of the structures table separated by tabs (\t)

structure_registry.tsv.gz            A compressed file of the structure_registry table separated by tabs (\t)

wurcs.tsv.gz                         A compressed file of the wurcs table separated by tabs (\t)

reference.tsv.gz                     A compressed file of the reference table separated by tabs (\t)

relation.tsv.gz                      A compressed file of the relation table separated by tabs (\t)

comments.tsv.gz                      A compressed file of the comments table separated by tabs (\t)

database_accession.tsv.gz            A compressed file of the database_accession table separated by tabs (\t)

names.tsv.gz                         A compressed file of the names table separated by tabs (\t)

source.tsv.gz                        A compressed file of the source table separated by tabs (\t)

status.tsv.gz                        A compressed file of the source table separated by tabs (\t)
