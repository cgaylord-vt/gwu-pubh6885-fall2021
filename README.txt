The dada2 to omeClust and MaAsLin2 pipeline applied to
cesarean birth seeding data.

Clark Gaylord
cgaylord@gwu.edu
December 2021

Example is to analyze one-month-old infants, fecal sample.
Other sites, times could also be done by querying the archive
for the accessions.

These are the files/directories you need to start:

-- ENAfiles 
   -- directory of ENA access files
   -- the search output and download script are provided

-- Metadata_Baby_Seeding_all_samples_final.xlsx 
   -- mess of metadata, use sheet labeled "Possible"
   -- this is being replaced with a text file fecal_onemonth-metadata.txt

-- RunningStart.R
   -- the script

-- we are removing the taxonomy component for the illustration
   -- silva_nr_v132_train_set.fa.gz
   -- silva_species_assignment_v132.fa.gz
   -- taxonomy data
