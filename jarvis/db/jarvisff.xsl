<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">

<html> 
<head>
<style>
.jcentered {
  margin: auto;
  width: 95%;
 border: 1px solid black;
  padding: 1px;
}
.jright {
  position: absolute;
  right: 0px;
  width: 300px;
  border: 3px solid #73AD21;
  padding: 10px;
}
.jleft {
  position: absolute;
  left: 0px;
  width: 300px;
  border: 3px solid #73AD21;
  padding: 10px;

}

</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.4.0/3Dmol-min.js"></script>
<script src="https://cdn.plot.ly/plotly-latest.min.js"> </script>
 


<!-- Bootstrap CSS START-->

</head>


<body>


  
<!-- Basic table start -->
  <table class="jcentered table  btn-table" >

 
    <tr>
      <td>ID: <xsl:value-of select="basic_info/id"/></td>
      <td>Pair style: <xsl:value-of select="basic_info/pair_style"/></td>
      <td>Input cell</td>
      <td>Input cell</td>
      <td>Optimized cell</td>
      <td>Optimized cell</td>
    </tr>
    <tr>
      <td>Chemical formula: <xsl:value-of select="basic_info/formula"/></td>
      <td>Pair-coeff.: <xsl:value-of select="basic_info/pair_coeff"/></td>
      <td>a <xsl:value-of select="basic_info/initial_a"/> &#8491;</td>
      <td>&#945;: <xsl:value-of select="basic_info/initial_alpha"/> &#176;</td>
      <td>a <xsl:value-of select="basic_info/final_a"/> &#8491;</td>
      <td>&#945;: <xsl:value-of select="basic_info/final_alpha"/> &#176;</td>
    </tr>
    <tr>
      <td>Input Space-group: <xsl:value-of select="basic_info/initial_spacegroup_symbol"/></td>
      <td>Final Space-group: <xsl:value-of select="basic_info/final_spacegroup_symbol"/></td>
      <td>b <xsl:value-of select="basic_info/initial_b"/> &#8491;</td>
      <td>&#946;: <xsl:value-of select="basic_info/initial_beta"/> &#176;</td>
      <td>b <xsl:value-of select="basic_info/final_b"/> &#8491;</td>
      <td>&#946;: <xsl:value-of select="basic_info/final_beta"/> &#176;</td>
    </tr>
    <tr>
      <td>Initial crystal system: <xsl:value-of select="basic_info/initial_crystal_system"/></td>
      <td>Final crystal system: <xsl:value-of select="basic_info/final_crystal_system"/></td>
      
      <td>c <xsl:value-of select="basic_info/initial_c"/> &#8491;</td>
      <td>&#947;: <xsl:value-of select="basic_info/initial_gamma"/> &#176;</td>
      <td>c <xsl:value-of select="basic_info/final_c"/> &#8491;</td>
      <td>&#947;: <xsl:value-of select="basic_info/final_gamma"/> &#176;</td>
    </tr>
    <tr>
      <td>Data source: <xsl:value-of select="basic_info/data_source"/></td>
      <td>Energy/atom (eV): <xsl:value-of select="basic_info/energy_per_atom"/></td>

      <td>Initial density (gcm<sup>-3</sup>): <xsl:value-of select="basic_info/initial_density"/></td>
      <td>Final density (<span>&#8491;</span><sup>3</sup>): <xsl:value-of select="basic_info/final_density"/></td>
      <td>Pressure: <xsl:value-of select="basic_info/pressure"/></td>
      <td>Number of species: <xsl:value-of select="basic_info/number_uniq_species"/></td>
    </tr>    
   
   
    




 
    
   
  </table>
  
  
<!-- Basic table end -->
<br></br>
<!-- Structure viewer start -->
<h3 style="text-align:center">Visualizing atomic structure</h3>
  <div >

  <div 
  class="jcentered"
  id="geometry"
  style="height: 400px; width: 400px; position: relative;"
  data-backgroundcolor="0xffffff"
  >  
  <span
    style="
      display: inline-block;
      z-index: 1;
      position: absolute;
      bottom: 0px;
      right: 0px;
      font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    "
    >
    
    </span>
  </div>

  
  

<script>
  var data =
    <xsl:value-of select="basic_info/xyz"/>;
  var viewer = $3Dmol.createViewer("geometry", {
    backgroundColor: "white",
  });
 // viewer.setStyle({stick:{colorscheme:"Jmol"}});
  
  viewer.addModel(data, "xyz");
  viewer.setStyle({ sphere: { radius: 0.7 } });
  
  viewer.addStyle({ stick: { radius: 0.2 } });
  viewer.zoomTo();
  viewer.setClickable({},true,function(atom,viewer,event,container) {
                       viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8});
                   });
 
  viewer.render();
</script>


  </div>
 <!-- Structure viewer end -->


 
 


  
  
  
  
  
  
  
  
  
  












<!--Elastic tensor table starts--> 


 
  <script>
  var x= <xsl:value-of select="basic_info/elastic_tensor/cij"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="jcentered";
   var text = document.createTextNode("Finite-difference Elastic constant-tensor");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
     var divElement = document.createElement("table");
    divElement.id = "elastic_tensor";
    divElement.setAttribute('class', 'jcentered');
    document.body.appendChild(divElement); 
    
    };
</script> 
<!--<table class="jcentered" style="border-spacing: 15px" id="elastic_tensor">  </table>
-->



<script> 
 var table = document.getElementById("elastic_tensor");
 elasticij=<xsl:value-of select="basic_info/elastic_tensor/cij"/> ;
 elasticij=elasticij.split(';');
 var arr=[["j1","j2","j3","j4","j5","j6"],];
 var row = table.insertRow(-1);
 var count=0;
 var tmp=[];
 for (var i=0;i&lt;elasticij.length;i++) {
 tmp=elasticij[i].split(',').map(Number);

 arr.push(tmp);
 

 };
   var  unit_system='GPa'; 
 var data = [{
  type: 'table',
 header:{values: [unit_system,"i1","i2","i3","i4","i5","i6"],height: 30,font: {family: "Arial", size: 24, color: ["black"]}},

  cells: {
    values: arr,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 30,
  }
}]






var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
}

Plotly.newPlot('elastic_tensor',data);
</script> 
<!--Elastic tensor table ends--> 
<br></br>

<!--Elastic tensor phdos and phbstructure starts--> 

 
<script>
  var x= <xsl:value-of select="basic_info/phonon_dos_line/phonon_dos_frequencies"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="jcentered";
   var text = document.createTextNode("FD Phonon DOS and Band structure");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    var divElement = document.createElement("div");
    divElement.id = "fddos";
    divElement.setAttribute('class', 'jcentered');
    document.body.appendChild(divElement);
    
    };
</script> 

 
<script>
var data=[];
function fddos_plotly(data){


var plot_font = 14;
var layout_convg = {
grid: {rows: 1, columns: 2},
  xaxis1: {domain: [0.1, 0.7],tickfont: {size: plot_font,color:'black'},title:{text: 'Frequency (cm-1)',font:{size: plot_font,color:"black"}}},
  yaxis1:{tickfont: {size: plot_font,color:'black'},title:{text: 'Density of states',font:{size: plot_font,color:"black"}}},
  yaxis2: {tickfont: {size: plot_font,color:'black'},anchor: 'x2',title:{text: 'Frequency (cm-1)',font:{size: plot_font,color:"black"}}},
  xaxis2: {tickfont: {size: plot_font,color:'black'},domain: [0.7, .99],title:{text: 'K-Point',font:{size: plot_font,color:"black"}}},
  
  showlegend: false,
  autosize: false,
     width:1400,
  
};

Plotly.newPlot('fddos', data, layout_convg);

};

  var energies= <xsl:value-of select="basic_info/phonon_dos_line/phonon_dos_frequencies"/>;
  energies=energies.split(',').map(Number);
  var ph_dos= <xsl:value-of select="basic_info/phonon_dos_line/phonon_dos_intensity"/>;
  ph_dos=ph_dos.split(',').map(Number);

  
  var data1 = {
    x: energies ,
    y:ph_dos,
    xaxis: "x1",
    yaxis: "y1",
    type: 'scatter',
    
  };

  data.push(data1);
  
   <!-- var x= <xsl:value-of select="basic_info/phonon_band_line/phonon_bandstructure_distances"/>; -->
  <!-- x=x.split(','); -->
  
  <!-- var y2= <xsl:value-of select="basic_info/phonon_band_line/phonon_bandstructure_frequencies"/>; -->
  <!-- y2=y2.split(';'); -->

  
  
<!-- var data_tmp = []; -->



    
 
  <!-- for (var i=0;i&lt;y2.length;i++) { -->
    
    <!-- var data2 = { -->
    <!-- x: x, -->
    <!-- y: y2[i].split(',').map(Number), -->
    <!-- marker:{color:'blue'}, -->
    <!-- showlegend: false, -->
    <!-- mode:'lines', -->
    <!-- xaxis: "x2", -->
    <!-- yaxis: "y2", -->
  
    <!-- type: 'scatter', -->
  <!-- }; -->
    <!-- data.push(data2); -->
  <!-- }; -->
  
    

  
  
  
  




  fddos_plotly(data);
</script>
<br></br>

<!--Elastic tensor phdos and phbstructure ends--> 




<!--Vacancy table starts--> 


 
  <script>
  var x= <xsl:value-of select="basic_info/vacancy_info"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="jcentered";
   var text = document.createTextNode("Vacancy formation energy");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
     var divElement = document.createElement("table");
    divElement.id = "vacancy";
    divElement.setAttribute('class', 'jcentered');
    document.body.appendChild(divElement); 
    
    };
</script> 


<script> 
 var table = document.getElementById("vacancy");
 var vac_data=<xsl:value-of select="basic_info/vacancy_info"/> ;
 vac_data=vac_data.split(';');
 var arr=[["Element","Wyck. mult.","Energy"],];
 var row = table.insertRow(-1);
 var count=0;
 var tmp=[];
 var header_arr=[]
  for (var i=0;i&lt;vac_data.length-1;i++) {
 tmp=vac_data[i].split(',');
   if (i===0){{header_arr.push('Atom') };
   }else{
   header_arr.push(i);
   };

 arr.push(tmp);
 

 };

  header_arr.push(vac_data.length-1);
 var data = [{
  type: 'table',
 header:{values: header_arr,height: 30,font: {family: "Arial", size: 24, color: ["black"]}},

  cells: {
    values: arr,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 24, color: ["black"]},
    height: 30,
  }
}]


var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
}

Plotly.newPlot('vacancy',data);
</script> 
<!--Vacancy table ends--> 
<br></br>






<!--Surface table starts--> 


 
  <script>
  var x= <xsl:value-of select="basic_info/surface_info"/>;
  if (x!==''){
   console.log('header test pass');
   var header = document.createElement("h3");
   header.className="jcentered";
   var text = document.createTextNode("Surface energy");
    header.appendChild(text);
    header.style.textAlign='center';
    document.body.appendChild(header);
    
    
     var divElement = document.createElement("table");
    divElement.id = "surface";
    divElement.setAttribute('class', 'jcentered');
    document.body.appendChild(divElement); 
    
    };
</script> 


<script> 
 var table = document.getElementById("vacancy");
 var vac_data=<xsl:value-of select="basic_info/surface_info"/> ;
 vac_data=vac_data.split(';');
 var arr=[["Miller","Energy"],];
 var row = table.insertRow(-1);
 var count=0;
 var tmp=[];
 var header_arr=[]
 
 for (var i=0;i&lt;vac_data.length-1;i++) {
 tmp=vac_data[i].split(',');
   if (i===0){{header_arr.push('Surface') };
   }else{
   header_arr.push(i);
   };

 arr.push(tmp);
 

 };
 header_arr.push(vac_data.length-1);
 var data = [{
  type: 'table',
 header:{values: header_arr,height: 30,font: {family: "Arial", size: 18, color: ["black"]}},

  cells: {
    values: arr,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 18, color: ["black"]},
    height: 30,
  }
}]


var layout = {
   xaxis1: {domain: [0.1, 0.5],},
   xaxis2: {domain: [0.5, 0.95],}, 
  
}

Plotly.newPlot('surface',data);
</script> 
<!--Surface table ends--> 
<br></br>





</body>
</html>
</xsl:template>
</xsl:stylesheet>



