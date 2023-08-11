document$.subscribe(function() {
  var tables = document.querySelectorAll("#j_table")
  tables.forEach(function(table) {
    new Tablesort(table)
  })
})

