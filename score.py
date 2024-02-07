import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from joblib import load
import argparse
from sklearn.metrics import precision_score, recall_score


def load_data(input_dir, sample_size):
    features_to_normalize = ['cov', 'abundance', 'count1', 'count2',
            'start_loss1', 'end_loss1', 'start_loss2', 'end_loss2', 
            'start_loss3', 'end_loss3','start_merged_loss', 'end_merged_loss',
            'seq_min_cnt','seq_min_abd', 'seq_max_cnt','seq_max_abd'
            ]    
    frames = []

    columns = [
        "tid", "meta_tid", "chr", "cov", "cov2", "abundance", "confidence",
        "count1", "count2", "num_exons", "gr_vertices", "gr_edges", "gr_reads",
        "gr_subgraph", "v", "e", "junc_ratio", "max_mid_exon_len", "start_loss1",
        "start_loss2", "start_loss3", "end_loss1", "end_loss2", "end_loss3",
        "start_merged_loss", "end_merged_loss", "introns", "intron_ratio",
        "start_introns", "start_intron_ratio", "end_introns", "end_intron_ratio",
        "uni_junc", "seq_min_wt", "seq_min_cnt", "seq_min_abd", "seq_min_ratio",
        "seq_max_wt", "seq_max_cnt", "seq_max_abd", "seq_max_ratio", "start_cnt",
        "start_weight", "start_abd", "end_cnt", "end_weight", "end_abd",
        "unbridge_start_coming_count", "unbridge_start_coming_ratio",
        "unbridge_end_leaving_count", "unbridge_end_leaving_ratio"
    ]

    # Load all data frames
    for i in range(0, sample_size+1):
        df = pd.read_csv(f"{input_dir}/{i}.trstFeature.csv", dtype={2: str}, names=columns, header=None, sep='\t')
        #df=pd.read_csv(f"{input_dir}/{i}.stats.csv", dtype={2: str})
        df['meta_only'] = (i == sample_size) and (df['count2'] == 1)
        df['sample_id'] = i
        df['sample_size'] = sample_size
        if(i == sample_size): df = df[df['count2'] == 1]
        frames.append(df)

    combined_df = pd.concat(frames, ignore_index=True)
    print(combined_df)

    max_cnt = combined_df['count2'].max()
    combined_df[features_to_normalize] = combined_df[features_to_normalize]/max_cnt
    print("max_cnt: ", max_cnt)
    print(combined_df)
    return combined_df


def main(args):
    print(f"Directory of individual samples: {args.input_dir}")
    print(f"Pre-trained model file: {args.model}")
    print(f"Number of samples: {args.count}")
    print(f"Output file: {args.output_file}")
    print(f"Min probability score: {args.prob_score}")

    model = load(args.model)

    # Test
    test_data = load_data(args.input_dir, args.count)

    X_test = test_data[['cov', 'cov2', 'abundance', 'confidence', 'count1', 'count2',
              'num_exons', 'gr_vertices', 'gr_edges', 'v', 'e', 
              'junc_ratio', 'max_mid_exon_len',
              'start_loss1', 'end_loss1', 'start_loss2', 'end_loss2', 
              'start_loss3', 'end_loss3','start_merged_loss', 'end_merged_loss',
              'introns', 'intron_ratio', 'start_introns', 'end_introns',
              'start_intron_ratio', 'end_intron_ratio','uni_junc', 
              'seq_min_wt','seq_min_cnt','seq_min_abd','seq_min_ratio', 
              'seq_max_wt','seq_max_cnt','seq_max_abd','seq_max_ratio', 
              'meta_only', 'sample_size',
              'start_cnt', 'start_weight', 'start_abd',
              'end_cnt', 'end_weight', 'end_abd',
              'gr_reads', 'gr_subgraph',
              "unbridge_start_coming_count", "unbridge_start_coming_ratio",
              "unbridge_end_leaving_count", "unbridge_end_leaving_ratio",
              ]]

    print(X_test.shape)

    y_prob = model.predict_proba(X_test)[:, 1]

    test_data['y_prob'] = y_prob
    meta_df = test_data.groupby('meta_tid').agg(
            y_prob_mean=('y_prob', 'mean'),
            #common_label=('Matched', 'first')  # All labels within a group are the same
            ).reset_index()

    output_df = meta_df[meta_df['y_prob_mean'] >= args.prob_score]
    print(output_df)
    output_df = output_df[['meta_tid', 'y_prob_mean']]
    print(output_df)
    output_df.to_csv(args.output_file, index=False)
    print(f"Save scores into {args.output_file}.")    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score input transcripts using a pre-trained model.")
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="Directory of results of individual samples")
    parser.add_argument("-m", "--model", type=str, required=True, help="Pre-trained model file")
    parser.add_argument("-c", "--count", type=int, required=True, help="Number of samples")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="Output file for scored transcripts")
    parser.add_argument("-p", "--prob_score", type=float, required=False, default=0.2, help="Minimum probability score(default = 0.2)")

    args = parser.parse_args()
    main(args)
