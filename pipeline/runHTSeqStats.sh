find . -name "*.txt" | xargs -I {} bash -c 'echo $0 $(grep __no_feature $0) $(grep __alignment_not_unique $0)' {}
find . -name "*.txt" -exec awk 'BEGIN{sum=0}{sum+=$2}END{print sum}' "{}" \;

