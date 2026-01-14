python CNN_Local_noAutoencoders_3x3.py \
  --working_dir "/your/working/directory" \
  --sdy569_feature "data/SDY569_tidy.csv" \
  --sdy797_feature "data/SDY797_tidy.csv" \
  --sdy1737_feature "data/SDY1737_tidy.csv" \
  --sdy569_cpeptide "data/SDY569_cpeptide_auc_tidy.csv" \
  --sdy797_cpeptide "data/SDY797_cpeptide_auc_tidy.csv" \
  --sdy1737_cpeptide "data/SDY1737_cpeptide_auc_tidy.csv" \
  --sdy569_cnn_test "test_data/SDY569_3x3_test.csv" \
  --sdy797_cnn_test "test_data/SDY797_3x3_test.csv" \
  --sdy1737_cnn_test "test_data/SDY1737_3x3_test.csv"

python CNN_Federated_noAutoencoder_3x3.py \
  --sdy569_train data/cleaned/SDY569_train.csv \
  --sdy797_train data/cleaned/SDY797_train.csv \
  --sdy1737_train data/cleaned/SDY1737_train.csv \
  --sdy569_test data/cleaned/SDY569_test.csv \
  --sdy797_test data/cleaned/SDY797_test.csv \
  --sdy1737_test data/cleaned/SDY1737_test.csv
