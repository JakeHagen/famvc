#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include <htslib/khash.h>
#include <htslib/vcfutils.h>
KHASH_SET_INIT_STR(strs)
KHASH_MAP_INIT_STR(strm, char*)

int parse_ped(char *fname, khash_t(strm) *h, khint_t *n_fids) {
    htsFile *fp = hts_open(fname, "r");
    //if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) exit(1);//error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL, ret;
    char *fid = NULL;
    char *sid = NULL;

    khash_t(strs) *fids;
    fids = kh_init(strs);
    khint_t k;
    
    do
    {
	int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<2 ) { break; }//error("Could not parse the ped file: %s\n", str.s);

        char *fid = &str.s[off[0]];
	char *sid = &str.s[off[1]];

	//put fid in set to later get unique fids
	kh_put(strs, fids, strdup(fid), &ret);

	//set sample as key, family as value
	k = kh_put(strm, h, strdup(sid), &ret);
	kh_val(h, k) = strdup(fid);

    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    *n_fids = kh_size(fids);
    kh_destroy(strs, fids);

    free(str.s);
    free(fid);
    free(sid);
    free(off);
    hts_close(fp);
    return 0;
}

int main(int argc, char *argv[]) {
    //create family sample map from pedigree

    khash_t(strm) *fsm = kh_init(strm);
    int ret;
    khint_t n_fids = 0;
    ret = parse_ped(argv[1], fsm, &n_fids);
    
    htsFile *fp = hts_open(argv[2], "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf_hdr_append(hdr,"##INFO=<ID=fam_VC,Number=1,Type=Integer,Description=\"\">");
    bcf_hdr_append(hdr,"##INFO=<ID=fam_VF,Number=1,Type=Float,Description=\"\">");

    bcf1_t *line = bcf_init();

    htsFile *out = hts_open("-", "wb");

    ret = bcf_hdr_write(out, hdr);
    if (ret!=0) {
	    printf("error writing header");
	    exit(ret);
    }

    khash_t(strs) *h;
    int absent;
    khint_t k;
    
    h = kh_init(strs);
    int ped_warn_c = 0;
    while (bcf_read(fp, hdr, line) == 0) {
	    bcf_fmt_t *gt = bcf_get_fmt(hdr, line, "GT");
	    int i;
	    for (i=0; i<line->n_sample; i++) {
		    int type = bcf_gt_type(gt,i,NULL,NULL);
		    if (type == 2 || type == 1) {
			    char *sid = hdr->samples[i];

			    k = kh_get(strm, fsm, sid);
			    //char *fid = kh_val(fsm, k);

			    
			    if (!kh_val(fsm, k)) {
				    if (ped_warn_c == 10) {
					fprintf( stderr, "warning %s not in pedigree, disregarding, WILL NOT PRINT MESSAGE ANYMORE\n", sid);
					ped_warn_c++;
				    }
				    if (ped_warn_c < 10) {
					    fprintf( stderr, "warning %s not in pedigree, disregarding\n", sid);
					    ped_warn_c++;
				    }
				    continue;
				    
			    }
			    kh_put(strs, h, kh_val(fsm, k), &absent);
		    }
	    }

	    khint_t famac = kh_size(h);
	    float famaf = (float) famac / (float) n_fids;
	    bcf_update_info_int32(hdr, line, "fam_VC", &famac, 1);
	    bcf_update_info_float(hdr, line, "fam_VF", &famaf, 1);

	    ret = bcf_write(out, hdr, line);
	    if (ret!=0) {
		    printf("error writing line");
		    exit(ret);
    	    }
	    kh_clear(strs, h);
    }

    for (k = 0; k < kh_end(fsm); ++k) {
	    if (kh_exist(fsm, k)) {
		    free((char*)kh_key(fsm, k));
		    free((char*)kh_val(fsm, k));
	    }
    }
    kh_destroy(strm, fsm);
    kh_destroy(strs, h);
    //khash_str2str_destroy_free_all(fam_samp);

    
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    hts_close(out);
    bcf_destroy(line);
}
