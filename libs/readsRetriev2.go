package main

import (
	"fmt"
	"os"
	"sort"
	"bufio"
	"flag"
	"strings"
	"strconv"

	"github.com/lskatz/FqReader.go"
)


func validBC(f string, cells int) map[string]int {
	fh, err := os.Open(f)
	if err != nil {
		fmt.Println("%v\n", err)
		os.Exit(1)
	}
	defer fh.Close()

	n := 0
	barcodes := make(map[string]int)
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() {
		if n > cells - 1 {
			break
		}

		arr := strings.Split(scanner.Text(), "\t")
		barcodes[arr[0]] = 1
		n += 1
	}
	fmt.Printf("barcode number: %v\n", len(barcodes))

	return barcodes
}

func main() {
	//start := time.Now()
	var fbc, cells, fq, prefix, umi_tab string
	flag.StringVar(&fbc, "bc", "" , "barcode list")
	flag.StringVar(&cells, "cells", "" , "cells number")
	flag.StringVar(&fq, "fq", "" , "fastq")
	flag.StringVar(&umi_tab, "tab", "" , "umi counts table")
	flag.StringVar(&prefix, "o", "" , "output prefix")
	flag.Parse()

	cellNum, err := strconv.Atoi(cells)
	var fqr FqReader.FqReader

	fout, _ := os.OpenFile(prefix + "_final.fq", os.O_CREATE|os.O_WRONLY, 0666)
	defer fout.Close()

	fhr, err := os.Open(fq)
	if err != nil {
		fmt.Println("%v\n", err)
		os.Exit(1)
	}
	defer fhr.Close()

	bc_keep := validBC(fbc, cellNum)
	umi_cts := make(map[string]int)

	fqr.Reader = bufio.NewReader(fhr)
	for r, done := fqr.Iter(); !done; r, done = fqr.Iter() {
		bc := r.Seq[0:20]
		if _, indb := bc_keep[bc]; indb {
                        // new_id := r.Name + "_" + bc_umi
			seq := r.Seq[28:]
			qual := r.Qual[28:]
			new_id := r.Name + "_" + bc
                        fout.WriteString(new_id+"\n"+seq+"\n+\n"+qual+"\n")
			umi_cts[bc] += 1
                } else {
                        continue
                }
	}

	if umi_tab == "on" {
		fout2, _ := os.OpenFile(prefix + "_UMI_counts.tsv", os.O_CREATE|os.O_WRONLY, 0666)
		defer fout2.Close()

		keys := make([]string, 0, len(umi_cts))

		for key := range umi_cts {
			keys = append(keys, key)
		}
		sort.SliceStable(keys, func(i, j int) bool{
			return umi_cts[keys[i]] > umi_cts[keys[j]]
		})

		for _, k := range keys{
			//fmt.Println(k, umi_cts[k])
			umi := strconv.Itoa(umi_cts[k])
			fout2.WriteString(k + "\t" + umi + "\n")
		}
	}
	//elapsed := time.Since(start)
	//fmt.Println("total time: ", elapsed)
}

