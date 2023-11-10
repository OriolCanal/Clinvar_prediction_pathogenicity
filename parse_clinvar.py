import pandas as pd
import gzip

def clinvar_vcf_to_pd(vcf_path):
    # list where we will store dictionary params of each variant
    variants_params = list()
    with gzip.open(vcf_path, "rt") as file:
        for line in file:
            # dictinary where we will store each parameter of the variant
            if line.startswith("#"):
                # descriptor lines not interested in
                continue
            fields = line.split("\t")
            # dictinary where we will store each parameter of the variant
            # obtaining parameters from each variant 
            chrom = fields[0]
            pos = fields[1]
            id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter = fields[6]
            info = fields[7]
            dict_params = {
                "Chrom" : chrom,
                "Pos" : pos,
                "Id" : id,
                "Ref" : ref,
                "Alt" : alt,
                "Qual" : qual,
                "Filter" : filter,

            }

            # in info we have different parameters 
            clnv_params = info.split(";")

            for clnv_param in clnv_params:
                key_value = clnv_param.split("=")
                key = key_value[0]
                value = key_value[1]
                # it's a comma seperated list of molecular consequences
                if key == "MC":
                    if "," in value:
                        mol_conseqs_ids = list()
                        mol_conseqs = value.split(",")
                        # print(mol_conseqs)
                        for mol_conseq in mol_conseqs:
                            # taking the string id to be converted into factors 
                            # print(mol_conseq)
                            mol_conseq_id = mol_conseq.split("|")[1]
                            mol_conseqs_ids.append(mol_conseq_id)
                        value = ",".join(mol_conseqs_ids)
                    else:
                        value = value.split("|")[1]


                dict_params[key] = value
            

            variants_params.append(dict_params)
    return(variants_params)


variants_params = clinvar_vcf_to_pd("/home/ocanal/ANN_DIR/clinvar/hg38/clinvar_20231104.vcf.gz")
df = pd.DataFrame(variants_params)




