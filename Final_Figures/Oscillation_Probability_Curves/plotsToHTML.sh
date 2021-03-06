#! /bin/sh


TNsize=600
nx=1

#for dir in `find ${destination} -maxdepth 1 -mindepth 1 -printf "%f\n"`
    
  echo "Directory ${dir}"

  mv pdf/* .

  rm -rf pdf small gif index.html

  mkdir -p pdf
  mv *.pdf pdf

  file="index.html"
  blue="#0000FF"
  red="#FF0000"
  black="#000000"
  title="${dir}"

  mkdir gif
  mkdir small

  rm -rf ${file}

  echo "<HTML>" >> ${file}
  echo "<TITLE>${title}</TITLE>" >> ${file}
  echo " <BODY>" >> ${file}
  # echo "  <CENTER><H1><font color=\"${blue}\">${title}</font></H1></CENTER>" >> ${file}
  echo "<p>" >> ${file}
  echo ${text} >> ${file}
  echo "</p>" >> ${file}

  echo "  <TABLE>" >> ${file}
  echo "   <TR>" >> ${file}

  # echo "<a href=\"Events.txt\"><H3><font color=\"${blue}\"> Events and scale </font></H3></a>" >> ${file}

  echo "Creating thumbnails..."

  #########################
  i1=0
  j1=0
  echo "    </TR>" >> ${file}
  echo "  </TABLE>" >> ${file}
  echo "<H1><font color=\"${blue}\">  </font></H1>" >> ${file}
  echo "  <TABLE>" >> ${file}
  echo "   <TR>" >> ${file}

  for pdffile1 in `cd pdf; ls *.pdf` ; do 

    echo "Converting pdf to gif, ${pdffile1} ..."

      pdf1=$pdffile1

      #creating scaled thumbnails !!!
      img1=`echo ${pdf1} | sed 's/.pdf/.gif/'`
      convert -density 300 -trim -quality 100 +repage pdf/${pdf1} gif/${img1}
      convert -scale ${TNsize} -trim gif/${img1} small/TN_${img1}

      echo "    <td align=\"center\">" >> ${file}
      echo "      <a href=\"gif/${img1}\"> <img src=\"small/TN_${img1}\"> </a>" >> ${file}
      pdf=`echo $img1 | sed 's/.gif/.pdf/'`
      echo "      <br><a href=\"gif/${img1}\">${img1}</a>, <a href=\"pdf/${pdf1}\">[pdf]</a>  "  >> ${file}
      echo "    </td>  "  >> ${file}
      i1=`expr $i1 + 1` 
      j1=`expr ${i1} % ${nx}`
      if [ $j1 -eq 0 ] ; then
         echo "   </TR><TR>" >> ${file}
      fi

  done   # loop over all files


  echo "    </TR>" >> ${file}
  echo "  </TABLE>" >> ${file}
  DATE=`date`
  echo "  <br><br>Created ${DATE} by (c) ${USER} using <a href=\"createHTMLpage.sh\">createHTMLpage.sh</a>"  >> ${file}
  echo " </BODY>" >> ${file}
  echo "</HTML>" >> ${file}


