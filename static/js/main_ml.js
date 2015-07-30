$(document).ready(function(){
       
    $('#txt_area').hide();
    $('#image').hide();
    
    $("#get_smiles").on('click',  function(){

      $('#txt_area').val(ChemDoodle.writeMOL(sketcher.getMolecule()));
      var smiles_mol  = $("#txt_area").val();
      $.ajax({
	    type: "GET",
	    url: "/ml_app/get_smiles",
	    data : {'smiles_mol': smiles_mol},
	    success: function(data)
	    {
	    $('#smiles').val(data.smiles);          
		
	    },
	    error: function(error) 
	    {
		console.log(error);
	    }
	});
	

   });

  $('#run').on('click', function(){
      var smiles = $('#smiles').val();
      var ml_options = $('#options').val();
      if (ml_options === 'mp'){
	  $('#warning').text('Please be patient. It may take several seconds to generate the result.');
	 } else{$('#warning').text('')};
      $.ajax({
	  type: "GET",
	  url: "/ml_app/predicted",
	  data: {'smiles': smiles, 'ml_options': ml_options},
	  success: function(data2){
	      $('#image').show();
	      $('#outcome').attr('src', data2.image);
	      $('#units').text(data2.units);	            
	  },
	  error: function(error){
	      console.log(error);
	  }
      });  
  });

})
