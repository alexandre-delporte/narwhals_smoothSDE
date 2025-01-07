for file in plot_fjords*.png; do
    new_name=$(echo "$file" | sed 's/plot_fjords\(.*\)\.png/plot_fjords_12h_12ID_4km\1.png/')
    mv "$file" "$new_name"
done

