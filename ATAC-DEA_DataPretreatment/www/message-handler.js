// This recieves messages of type "testmessage" from the server.
// See http://shiny.rstudio.com/gallery/server-to-client-custom-messages.html
// for details
Shiny.addCustomMessageHandler("Message",
  function(message) {
    alert(JSON.stringify(message));
  }
);
