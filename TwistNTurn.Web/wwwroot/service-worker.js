// TwistNTurn Service Worker
const CACHE_NAME = 'twistnturn-v5';

// Install event: remain waiting until user accepts update via SKIP_WAITING
self.addEventListener('install', event => {
    // Do not call self.skipWaiting() unconditionally so active sessions keep matching assets
});

// Activate event: clean up outdated caches and take control
self.addEventListener('activate', event => {
    event.waitUntil(
        caches.keys().then(keys => Promise.all(
            keys.filter(k => k !== CACHE_NAME).map(k => caches.delete(k))
        )).then(() => self.clients.claim())
    );
});

// Fetch event:
// - Network-first for navigation, HTML, un-fingerprinted CSS, and manifest
// - Cache-first with offline fallback for fingerprinted _framework and image assets
self.addEventListener('fetch', event => {
    if (event.request.method !== 'GET') return;

    const url = new URL(event.request.url);

    // Only cache requests from our origin
    if (url.origin !== self.location.origin) return;

    const isNetworkFirst =
        event.request.mode === 'navigate' ||
        url.pathname.endsWith('/') ||
        url.pathname.endsWith('.html') ||
        url.pathname.endsWith('.css') ||
        url.pathname.endsWith('.webmanifest');

    if (isNetworkFirst) {
        event.respondWith(
            fetch(event.request)
                .then(networkResponse => {
                    if (networkResponse && networkResponse.status === 200) {
                        const copy = networkResponse.clone();
                        event.waitUntil(caches.open(CACHE_NAME).then(cache => cache.put(event.request, copy)).catch(() => {}));
                    }
                    return networkResponse;
                })
                .catch(() => caches.match(event.request).then(cached => cached || caches.match('./')).then(res => res || new Response('', { status: 503 })))
        );
        return;
    }

    // Fingerprinted / static assets (_framework, images): Cache-first with graceful offline catch
    event.respondWith(
        caches.match(event.request).then(cachedResponse => {
            if (cachedResponse) {
                return cachedResponse;
            }
            return fetch(event.request)
                .then(networkResponse => {
                    if (networkResponse && networkResponse.status === 200) {
                        const copy = networkResponse.clone();
                        event.waitUntil(caches.open(CACHE_NAME).then(cache => cache.put(event.request, copy)).catch(() => {}));
                    }
                    return networkResponse;
                })
                .catch(() => caches.match(event.request).then(cached => cached || new Response('', { status: 503 })));
        })
    );
});

// Message listener for skip waiting or cache busting
self.addEventListener('message', event => {
    if (event.data) {
        if (event.data.type === 'SKIP_WAITING') {
            self.skipWaiting();
        } else if (event.data.type === 'CLEAR_CACHE') {
            event.waitUntil(
                caches.keys().then(keys => Promise.all(keys.map(k => caches.delete(k))))
            );
        }
    }
});
