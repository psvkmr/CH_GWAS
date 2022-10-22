# AutoConvert using cmd windows
"C:\Program Files\Illumina\AutoConvert 2.0\AutoConvert.exe" I:\psivakumar\Emer_CH\Plate_3_IDATs\203692670081 I:\psivakumar\Emer_CH\Plate_3_GTC I:\psivakumar\GSA_v2_files\GSA-24v2-0_A1.bpm I:\psivakumar\GSA_v2_files\GSA-24v2-0_A1_ClusterFile.egt


python3.6 /data/kronos/Genetics_Software/GTCtoVCF/gtc_to_vcf.py \
--gtc-paths /array/psivakumar/Emer_CH/Plate_3_GTC/ \
--manifest-file /array/psivakumar/GSA_v2_files/GSA-24v2-0_A1.csv \
--genome-fasta-file /data/kronos/NGS_Reference/fasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
--output-vcf-path /array/psivakumar/Emer_CH/Plate_3_VCF/ \
--log-file /array/psivakumar/Emer_CH/Plate_3_VCF/GtcToVcf.log

# need to compress vcf with bgzip
bcftools merge 203692670068_R01C01.vcf 203692670068_R03C01.vcf -O v -o test.vcf
