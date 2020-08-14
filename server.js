const express = require("express");
const session = require("express-session");
const crypto = require("crypto");
const SpotifyWebApi = require("spotify-web-api-node");
const app = express();
const fetch = require("node-fetch");
var bodyParser = require("body-parser");
app.use(bodyParser.json()); // to support JSON-encoded bodies
app.use(
  bodyParser.urlencoded({
    // to support URL-encoded bodies
    extended: true
  })
);
app.use("/static", express.static("./static/"));

app.set("view engine", "pug");
app.use(express.static("public"));
app.set("trust proxy", 1);
app.use(
  session({
    resave: true,
    saveUninitialized: false,
    secret: process.env.SESSION_SECRET,
    maxAge: 60000 * 5
  })
);
app.get("/", function(request, response) {
  response.render("index");
});

const redirectUri = "https://track-visualizer.glitch.me/callback";
const scopes = [];
const showDialog = true;

let spotifyApi = new SpotifyWebApi({
  clientId: process.env.CLIENT_ID,
  clientSecret: process.env.CLIENT_SECRET,
  redirectUri: redirectUri
});

app.get("/authorize", function(request, response) {
  let state = crypto.randomBytes(12).toString("hex");
  request.session.state = state;
  let authorizeURL = spotifyApi.createAuthorizeURL(scopes, state, showDialog);
  response.redirect(authorizeURL);
});

app.get("/callback", function(request, response) {
  if (request.session.state !== request.query.state) {
    response.sendStatus(401);
  }
  let authorizationCode = request.query.code;
  spotifyApi.authorizationCodeGrant(authorizationCode).then(
    data => {
      request.session.access_token = data.body["access_token"];
      response.redirect("/playlists");
    },
    error => {
      console.log(
        "Something went wrong when retrieving the access token!",
        error.message
      );
    }
  );
});

app.get("/logout", function(request, response) {
  request.session.destroy();
  response.redirect("/");
});

app.get("/playlists", function(request, response) {
  let loggedInSpotifyApi = new SpotifyWebApi();
  loggedInSpotifyApi.setAccessToken(request.session.access_token);
  loggedInSpotifyApi.getUserPlaylists().then(
    data => {
      response.render("playlists", { playlists: data.body.items });
    },
    error => {
      response.sendStatus(400);
    }
  );
});

app.get("/playlists/:id", function(req, res) {
  let loggedInSpotifyApi = new SpotifyWebApi();
  loggedInSpotifyApi.setAccessToken(req.session.access_token);
  loggedInSpotifyApi.getPlaylist(req.params.id).then(
    data => {
      res.render("playlist", {
        tracks: data.body.tracks.items,
        name: data.body.name
      });
    },
    error => {
      res.sendStatus(400);
    }
  );
});

app.post("/playlists/:id", function(req, res) {
  const artist = replaceSpaces(req.body.artist_name);
  const song_name = replaceSpaces(req.body.song_title);
  console.log(
    "https://api.musixmatch.com/ws/1.1/matcher.lyrics.get?q_track=" +
      song_name +
      "&q_artist=" +
      artist +
      "&apikey=" +
      process.env.MUSIXMATCH_API_SECRET
  );
  fetch(
    "https://api.musixmatch.com/ws/1.1/matcher.lyrics.get?q_track=" +
      song_name +
      "&q_artist=" +
      artist +
      "&apikey=" +
      process.env.MUSIXMATCH_API_SECRET
  ).then(
    data => {
      data.json().then(parsed => {
        let trimmed_lyrics;
        if (parsed.message.body.lyrics === undefined) {
          trimmed_lyrics = "";
        } else {
          trimmed_lyrics = parsed.message.body.lyrics.lyrics_body.slice(0, -75);
        }
        res.render("lyrics", {
          song_title: req.body.song_title,
          track_id: req.body.track_id,
          playlist_id: req.params.id,
          lyrics: trimmed_lyrics
        });
      });
    },
    error => {
      console.log(error);
    }
  );
});

function replaceSpaces(string) {
  return string.replace(/\s/g, "%20");
}

let listener = app.listen(process.env.PORT, function() {
  console.log("Your app is listening on port " + listener.address().port);
});
