//#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "khash_str2str.h"
#include <htslib/kseq.h>
#include <htslib/khash.h>
#include <htslib/vcfutils.h>
KHASH_SET_INIT_STR(str)

int parse_ped(char *fname, void *hash)
{
    htsFile *fp = hts_open(fname, "r");
    //if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) exit(1);//error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    char *family = NULL;
    char *sample = NULL;
    do
    {
	int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<2 ) { break; }//error("Could not parse the ped file: %s\n", str.s);

        char *family = &str.s[off[0]];
	char *sample = &str.s[off[1]];

	//set sample as key, family as value
	khash_str2str_set(hash,strdup(sample),strdup(family));

    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );


    free(str.s);
    free(family);
    free(sample);
    free(off);
    hts_close(fp);
    return 0;
}

int main(int argc, char *argv[]) {
    void *fam_samp = khash_str2str_init();
    parse_ped(argv[1], fam_samp);
    htsFile *fp = hts_open(argv[2], "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf_hdr_append(hdr,"##INFO=<ID=fam_AC,Number=1,Type=Integer,Description=\"\">");

    bcf1_t *line = bcf_init();

    htsFile *out = hts_open("-", "wb");

    int ret;
    ret = bcf_hdr_write(out, hdr);
    if (ret!=0) {
	    printf("error writing header");
	    exit(ret);
    }

    khash_t(str) *h;
    int absent;
    h = kh_init(str);
    int ped_warn_c = 0;
    while (bcf_read(fp, hdr, line) == 0) {
	    bcf_fmt_t *gt = bcf_get_fmt(hdr, line, "GT");
	    int i;
	    int c = 0;
	    for (i=1; i<line->n_sample; i++) {
		    int type = bcf_gt_type(gt,i,NULL,NULL);
		    if (type == 2 || type == 1) {
			    //khash_str2str_set(hash,strdup(sample),strdup(family));
			    char *samp = hdr->samples[i];
			    char *fam = khash_str2str_get(fam_samp, samp);
			    if (!fam) {
				    if (ped_warn_c == 10) {
					fprintf( stderr, "warning %s not in pedigree, disregarding, WILL NOT PRINT MESSAGE ANYMORE\n", samp);
					ped_warn_c++;
				    }
				    if (ped_warn_c < 10) {
					    fprintf( stderr, "warning %s not in pedigree, disregarding\n", samp);
					    ped_warn_c++;
				    }
				    continue;
			    }
			    c++;
			    kh_put(str, h, fam, &absent);
		    }
	    }
	    int32_t famac = kh_size(h);
	    bcf_update_info_int32(hdr, line, "fam_AC", &famac, 1);
	    ret = bcf_write(out, hdr, line);
	    if (ret!=0) {
		    printf("error writing line");
		    exit(ret);
    	    }
	    kh_clear(str, h);
    }

    kh_destroy(str, h);
    khash_str2str_destroy_free_all(fam_samp);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    hts_close(out);
    bcf_destroy(line);
}
